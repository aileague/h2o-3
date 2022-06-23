package hex.tree.sdt;

import hex.ModelBuilder;
import hex.ModelCategory;
import hex.tree.sdt.binning.BinAccumulatedStatistics;
import hex.tree.sdt.binning.BinningStrategy;
import hex.tree.sdt.binning.Histogram;
import hex.tree.sdt.mrtasks.CountSplitValuesMRTask;
import hex.tree.sdt.mrtasks.GetClassCountsMRTask;
import hex.tree.sdt.mrtasks.SplitFrameMRTask;
import org.apache.log4j.Logger;
import water.DKV;
import water.exceptions.H2OModelBuilderIllegalArgumentException;
import water.fvec.Frame;
import water.util.*;

import java.lang.reflect.MalformedParametersException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Single Decision Tree
 */
public class SDT extends ModelBuilder<SDTModel, SDTModel.SDTParameters, SDTModel.SDTOutput> {
    private int _maxDepth;
    int _nodesCount;

    private double[][] _tree;
//    private double[][] _iterativelyBuiltTree;

    private Integer _actualDepth;

    /**
     * Tree root - used only when tree is built recursively
     */
    private Node _root;

    private SDTModel _model;
    transient Random _rand;

    // todo - create file with constants ?
    private final static int LIMIT_NUM_ROWS_FOR_SPLIT = 2; // todo - make a parameter with default value
    private final static double EPSILON = 0.000001d; 
    

    private static final Logger LOG = Logger.getLogger(SDT.class);


    public SDT(SDTModel.SDTParameters parameters) {
        super(parameters);
        _maxDepth = parameters._maxDepth;
        _actualDepth = 0;
        _nodesCount = 0;
        _tree = null;
        init(false);
    }

    public SDT(boolean startup_once) {
        super(new SDTModel.SDTParameters(), startup_once);
    }

    public double[][] compress() {
        //  if parent node is at index i in the array then the left child of that node is at index (2*i + 1) 
        //  and right child is at index (2*i + 2) in the array. 
        System.out.println("Nodes count when compressing: " + _nodesCount);
        
        // When built iteratively, tree is already in compressed array
        if(_tree != null) {
            return _tree;
        }
        // save linked tree to array
        // 2^k - 1 is max count of nodes, where k is depth
        _tree = new double[(int) Math.pow(2, _maxDepth + 1)][2];
        writeSubtreeStartingFromIndex(_root, 0);
        return _tree;
    }

    private void writeSubtreeStartingFromIndex(final Node actualNode, final int actualIndex) {
        if (actualNode == null) {
            return;
        }
        _tree[actualIndex][0] = actualNode.getFeature() == null ? -1
                : actualNode.getFeature().doubleValue();
        _tree[actualIndex][1] = actualNode.getThreshold() == null ? actualNode.getDecisionValue()
                : actualNode.getThreshold();
        writeSubtreeStartingFromIndex(actualNode.getLeft(), 2 * actualIndex + 1);
        writeSubtreeStartingFromIndex(actualNode.getRight(), 2 * actualIndex + 2);
    }

    public double[][] getCompressedTree() {
        if (_tree == null) {
            compress();
        }
        return _tree;
    }

    /**
     * This method is used only when we split data in each node (old implementation of findBestSplit)
     *
     * @param data - frame with relevant data for current node
     * @param featuresLimits - limits for each feature - in this case used for optimalization only
     * @return split info - holds feature index and threshold
     */
    private SplitInfo findBestSplit(final Frame data, DataFeaturesLimits featuresLimits) {
        // find split (feature and threshold)
        int featuresNumber = data.numCols();
        Pair<Double, Double> currentMinEntropyPair = new Pair<>(-1., Double.MAX_VALUE);
        int bestFeatureIndex = -1;
        for (int featureIndex = 0; featureIndex < featuresNumber - 1 /*last column is prediction*/; featureIndex++) {
            final int featureIndexForLambda = featureIndex;
            // iterate all candidate values of threshold
            Pair<Double, Double> minEntropyForFeature = featuresLimits.getFeatureRange(featureIndex)
                    .map(candidateValue -> new Pair<>(
                            candidateValue, calculateEntropyOfSplit(data, featureIndexForLambda, candidateValue)))
                    .min(Comparator.comparing(Pair::_2))
                    .get();
            if (minEntropyForFeature._2() < currentMinEntropyPair._2()) {
                currentMinEntropyPair = minEntropyForFeature;
                bestFeatureIndex = featureIndex;
            }
        }
        double threshold = currentMinEntropyPair._1();
        return new SplitInfo(bestFeatureIndex, threshold);
    }

    /**
     * When we don't split the frame, but use binning and update features limits for child nodes.
     *
     * @param histogram - histogram for relevant data
     * @return split info - holds feature index and threshold, null if the split could not be found.
     */
    private SplitInfo findBestSplit(Histogram histogram) {
        int featuresNumber = histogram.featuresCount();
        Pair<Double, Double> currentMinEntropyPair = new Pair<>(-1., Double.MAX_VALUE);
        int bestFeatureIndex = -1;
        for (int featureIndex = 0; featureIndex < featuresNumber - 1 /*last column is prediction*/; featureIndex++) {
            // skip constant features
            if (histogram.isConstant(featureIndex)) {
                continue;
            }
            // iterate all bins
            Pair<Double, Double> minEntropyForFeature = histogram.calculateBinsStatisticsForFeature(featureIndex).stream()
                    .filter(binStatistics -> ((binStatistics._leftCount >= LIMIT_NUM_ROWS_FOR_SPLIT)
                            && (binStatistics._rightCount >= LIMIT_NUM_ROWS_FOR_SPLIT))) // todo - consider setting min count of samples in bin 
//                    .peek(binStatistics -> {
//                        System.out.println("counts: " + binStatistics._maxBinValue + " " + binStatistics._leftCount + " " + binStatistics._rightCount);
//                    })
                    .map(binStatistics -> new Pair<>(
                            binStatistics._maxBinValue, calculateEntropyOfSplit(binStatistics)))
                    .min(Comparator.comparing(Pair::_2))
                    .orElse(null);
            if(minEntropyForFeature == null) {
                continue; // no split was found for this feature
            }
            if (minEntropyForFeature._2() < currentMinEntropyPair._2()) {
                currentMinEntropyPair = minEntropyForFeature;
                bestFeatureIndex = featureIndex;
            }
        }
        if(bestFeatureIndex == -1) {
            return null; // no split could be found
        }
        double threshold = currentMinEntropyPair._1();
        return new SplitInfo(bestFeatureIndex, threshold);
    }

    private Double binaryEntropy(int leftCount, int leftCount0, int rightCount, int rightCount0) {
        double a1 = (entropyBinarySplit(leftCount0 * 1.0 / leftCount)
                * leftCount / (leftCount + rightCount));
        double a2 = (entropyBinarySplit(rightCount0 * 1.0 / rightCount)
                * rightCount / (leftCount + rightCount));
        double value = a1 + a2;
//        System.out.println("value: " + value + ", t.l " + task.countLeft + ", t.l0 " + task.countLeft0 + ", t.r " + task.countRight + ", t.r0 " + task.countRight0);
        return value;
    }

    private double entropyBinarySplit(final double oneClassFrequency) { 
        return -1 * ((oneClassFrequency < EPSILON ? 0 : (oneClassFrequency * Math.log(oneClassFrequency)))
                + ((1 - oneClassFrequency) < EPSILON ? 0 : ((1 - oneClassFrequency) * Math.log(1 - oneClassFrequency))));
    }

    private Double calculateEntropyOfSplit(BinAccumulatedStatistics binStatistics) {
        return binaryEntropy(binStatistics._leftCount, binStatistics._leftCount0,
                binStatistics._rightCount, binStatistics._rightCount0);
    }

    public double calculateEntropyOfSplit(final Frame data, final int featureIndex, final double threshold) {
        CountSplitValuesMRTask task = new CountSplitValuesMRTask(featureIndex, threshold);
        task.doAll(data);
        return binaryEntropy(task.countLeft, task.countLeft0, task.countRight, task.countRight0);
    }

    /**
     * Select decision value for list.
     * @param zeroRatio #zero/#one in list
     * @return decision value (in current case - 0 or 1)
     */
    private int selectDecisionValue(double zeroRatio) {
        if (zeroRatio >= 0.5) {
            return 0;
        } else {
            return 1;
        }
    }

    /**
     * When tree is build recursively.
     * @param zeroRatio #zero/#one in list
     * @param node node to set decision value
     * @return list node
     */
    public Node makeListFromNode(double zeroRatio, Node node) {
        node.setDecisionValue(selectDecisionValue(zeroRatio));
        return node;
    }

    /**
     * When tree is build iteratively.
     * @param zeroRatio #zero/#one in list
     * @param nodeIndex node index
     */
    public void makeListFromNode(double zeroRatio, int nodeIndex) {
        _tree[nodeIndex][0] = -1; // list
        _tree[nodeIndex][1] = selectDecisionValue(zeroRatio);
        // nothing to return, node is modified inplace
    }

    /**
     * Log with base 2.
     * @param value value
     * @return log_2(value)
     */
    private int log2(int value) {
        return (int)(Math.log(value) / Math.log(2));
    }

    /**
     * Build ext node when tree is built iteratively. The queue is updated here.
     * @param limitsQueue queue with feature limits for nodes
     * @param nodeIndex index of node in the tree array
     */
    public void buildNextNode(Queue<DataFeaturesLimits> limitsQueue, int nodeIndex) {
        // take limits for actual node
        DataFeaturesLimits actualLimits = limitsQueue.poll();
        // if the element is null, then the node should not be built. Nulls exist to keep the array building straightforward
        if(actualLimits == null) {
            // don't save anything to tree (no node is created)
            // add imaginary left and right children to imitate right tree structure
            // left child
            limitsQueue.add(null);
            // right child
            limitsQueue.add(null);
            return;
        }
        
        // todo - add limit by information gain (at least because of ideal split for example 11111)
        // (count0, count1)
        Pair<Integer, Integer> classesCount = countClasses(actualLimits);
        double zeroRatio = classesCount._1() /*0*/ * 1.0 / (classesCount._1() /*0*/ + classesCount._2() /*1*/);
        if(nodeIndex == 1) {
            System.out.println("Classes counts in dataset: 0 - " + classesCount._1() +", 1 - " + classesCount._2());
        }
        // compute node depth
        int nodeDepth = (int) Math.floor(log2(nodeIndex + 1));
        if ((nodeDepth >= _maxDepth) || (classesCount._1() <= LIMIT_NUM_ROWS_FOR_SPLIT) || (classesCount._2() <= LIMIT_NUM_ROWS_FOR_SPLIT)
//                || zeroRatio > 0.999 || zeroRatio < 0.001
        ) {
//            System.out.println("Reason: depth=" + _actualDepth + ", ratio=" + zeroRatio + ", class0=" + classesCount._1() + ", class1=" + classesCount._2());
            
            // add imaginary left and right children to imitate valid tree structure
            // left child
            limitsQueue.add(null);
            // right child
            limitsQueue.add(null);
            makeListFromNode(zeroRatio, nodeIndex);
            return;
        }
        Histogram histogram = new Histogram(_train, actualLimits, BinningStrategy.EQUAL_WIDTH);

        SplitInfo bestSplitInfo = findBestSplit(histogram);
        // if no split could be found, make a list from current node
        if(bestSplitInfo == null) {
            // add imaginary left and right children to imitate right tree structure
            // left child
            limitsQueue.add(null);
            // right child
            limitsQueue.add(null);
            makeListFromNode(zeroRatio, nodeIndex);
            return;
        }

        _tree[nodeIndex][0] = bestSplitInfo._splitFeatureIndex;
        _tree[nodeIndex][1] = bestSplitInfo._threshold;

        DataFeaturesLimits limitsLeft = actualLimits.updateMax(bestSplitInfo._splitFeatureIndex, bestSplitInfo._threshold);
        // set new min to something bigger than threshold as the threshold is included in the left split and excluded in the right
        DataFeaturesLimits limitsRight = actualLimits.updateMin(bestSplitInfo._splitFeatureIndex, bestSplitInfo._threshold/* + 0.0000000001*/);
//        System.out.println("root: " + countClasses(actualLimits) + ", left: " + countClasses(limitsLeft) +
//                ", right: " + countClasses(limitsRight) + ", best feature: " + bestSplitInfo._splitFeatureIndex
//                + ", threshold: " + bestSplitInfo._threshold);

//        System.out.println("feature: " + bestSplitInfo._splitFeatureIndex + ", threshold: " + bestSplitInfo._threshold);
//        System.out.println("Left min-max: " + limitsLeft.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._min +
//                " " + limitsLeft.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._max);
//        System.out.println("Right min-max: " + limitsRight.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._min +
//                " " + limitsRight.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._max);
        
        // store limits for left child
        limitsQueue.add(limitsLeft);
        // store limits for right child
        limitsQueue.add(limitsRight);
    }


    /**
     * Used only when we split dataframes (old implementation of buildSubtree)
     * @param data
     * @param featuresLimits
     * @param nodeDepth
     * @return
     */
    public Node buildSubtree(final Frame data, DataFeaturesLimits featuresLimits, int nodeDepth) {
        Node subtreeRoot = new Node();
        _nodesCount++;
        if (_actualDepth < nodeDepth) {
            _actualDepth = nodeDepth;
        }
        // todo - add limit by information gain (at least because of ideal split for example 11111)
        double zeroRatio = getZeroRatio(data);
        if (nodeDepth >= _maxDepth || data.numRows() <= LIMIT_NUM_ROWS_FOR_SPLIT || zeroRatio > 0.9 || zeroRatio < 0.1) {
//            System.out.println("actualDepth: " + actualDepth + ", data.numRows(): " + data.numRows() + ", zeroRatio: " + zeroRatio);
            if (zeroRatio >= 0.5) {
                subtreeRoot.setDecisionValue(0);
            } else if (zeroRatio < 0.5) {
                subtreeRoot.setDecisionValue(1);
            }
            return subtreeRoot;
        }

        SplitInfo bestSplitInfo = findBestSplit(data, featuresLimits);


        subtreeRoot.setFeature(bestSplitInfo._splitFeatureIndex);
        subtreeRoot.setThreshold(bestSplitInfo._threshold);

        // split data
        DataSplit split = splitData(bestSplitInfo._splitFeatureIndex, bestSplitInfo._threshold);
//        String[] outputColNamesLeft = Arrays.stream(_train.names()).map(n -> n + "Left").toArray(String[]::new);
//        String[] outputColNamesRight = Arrays.stream(_train.names()).map(n -> n + "Right").toArray(String[]::new);
        subtreeRoot.setLeft(buildSubtree(split.leftSplit,
                featuresLimits.updateMax(bestSplitInfo._splitFeatureIndex, bestSplitInfo._threshold), nodeDepth + 1));
        subtreeRoot.setRight(buildSubtree(split.rightSplit,
                featuresLimits.updateMin(bestSplitInfo._splitFeatureIndex, bestSplitInfo._threshold), nodeDepth + 1));
        return subtreeRoot;
    }

    /**
     * Used when we use binning, but the tree is built recursively
     * @param featuresLimits
     * @param nodeDepth
     * @return
     */
    public Node buildSubtree(DataFeaturesLimits featuresLimits, int nodeDepth) {
        Node subtreeRoot = new Node();
        _nodesCount++;
        if (_actualDepth < nodeDepth) {
            _actualDepth = nodeDepth;
        }
        // todo - add limit by information gain (at least because of ideal split for example 11111)
        // (count0, count1)
        Pair<Integer, Integer> classesCount = countClasses(featuresLimits);
        double zeroRatio = classesCount._1() /*0*/ * 1.0 / (classesCount._1() /*0*/ + classesCount._2() /*1*/);
        if(nodeDepth == 1) {
            System.out.println("Classes counts in dataset: 0 - " + classesCount._1() +", 1 - " + classesCount._2());
        }
        if ((nodeDepth >= _maxDepth) || (classesCount._1() <= LIMIT_NUM_ROWS_FOR_SPLIT) || (classesCount._2() <= LIMIT_NUM_ROWS_FOR_SPLIT)
//                || zeroRatio > 0.999 || zeroRatio < 0.001
        ) {
//            System.out.println("Reason: depth=" + _actualDepth + ", ratio=" + zeroRatio + ", class0=" + classesCount._1() + ", class1=" + classesCount._2());
            return makeListFromNode(zeroRatio, subtreeRoot);
        }
        Histogram histogram = new Histogram(_train, featuresLimits, BinningStrategy.EQUAL_WIDTH);

        SplitInfo bestSplitInfo = findBestSplit(histogram);
        // if no split could be found, make a list from current node
        if(bestSplitInfo == null) {
            return makeListFromNode(zeroRatio, subtreeRoot);
        }

        subtreeRoot.setFeature(bestSplitInfo._splitFeatureIndex);
        subtreeRoot.setThreshold(bestSplitInfo._threshold);

        DataFeaturesLimits limitsLeft = featuresLimits.updateMax(bestSplitInfo._splitFeatureIndex, bestSplitInfo._threshold);
        // set new min to something bigger than threshold as the threshold is included in the left split and excluded in the right
        DataFeaturesLimits limitsRight = featuresLimits.updateMin(bestSplitInfo._splitFeatureIndex, bestSplitInfo._threshold/* + 0.0000000001*/);
//        System.out.println("root: " + countClasses(featuresLimits) + ", left: " + countClasses(limitsLeft) +
//                ", right: " + countClasses(limitsRight) + ", best feature: " + bestSplitInfo._splitFeatureIndex 
//                + ", threshold: " + bestSplitInfo._threshold);
        
//        System.out.println("feature: " + bestSplitInfo._splitFeatureIndex + ", threshold: " + bestSplitInfo._threshold);
//        System.out.println("Left min-max: " + limitsLeft.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._min +
//                " " + limitsLeft.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._max);
//        System.out.println("Right min-max: " + limitsRight.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._min +
//                " " + limitsRight.getFeatureLimits(bestSplitInfo._splitFeatureIndex)._max);
        subtreeRoot.setLeft(buildSubtree(limitsLeft, nodeDepth + 1));
        subtreeRoot.setRight(buildSubtree(limitsRight, nodeDepth + 1));

        return subtreeRoot;
    }


    /**
     * Compute initial features limits.
     * @return
     */
    private DataFeaturesLimits getInitialFeaturesLimits() {
        return new DataFeaturesLimits(
                Arrays.stream(_train.vecs())
                        // decrease min as the minimum border is always excluded and real min value could be lost
                        .map(v -> new FeatureLimits(v.min() - EPSILON, v.max()))
                        .collect(Collectors.toList()));
    }


    private class SDTDriver extends Driver {

        @Override
        public void computeImpl() {
            _model = null;
            try {
                init(false);
                if (error_count() > 0) {
                    throw H2OModelBuilderIllegalArgumentException.makeFromBuilder(SDT.this);
                }
                _rand = RandomUtils.getRNG(_parms._seed);
//                isolationTreeStats = new IsolationTreeStats();
                _model = new SDTModel(dest(), _parms,
                        new SDTModel.SDTOutput(SDT.this));
                _model.delete_and_lock(_job);
                buildSDT();
//                model._output._model_summary = createModelSummaryTable();
                LOG.info(_model.toString());
            } finally {
                if (_model != null)
                    _model.unlock(_job);
            }
        }

        /**
         * Select strategy for building sdt - splitting frame, recursively binning or iterative binning
         */
        private void buildSDT() {
            // switch to using or not using binning
//            if(Objects.equals(_parms._buildingStrategy, "splitting")) {
//                _root = buildSubtree(_train, getInitialFeaturesLimits(), 1); // without histogram
//            } else if(Objects.equals(_parms._buildingStrategy, "binning")) {
//                _root = buildSubtree(getInitialFeaturesLimits(), 0); // with histogram
//            } else if(Objects.equals(_parms._buildingStrategy, "iterative")){
//                buildSDTIteratively();
//            } else {
//                System.out.println("No building strategy found for name " + _parms._buildingStrategy);
//                throw new MalformedParametersException();
//            }
//            root = buildSubtree(_train, getInitialFeaturesLimits(), 1); // without histogram
//            _root = buildSubtree(getInitialFeaturesLimits(), 0); // with histogram // todo - 0 or 1 (0 works fine as well)
            buildSDTIteratively();
            System.out.println("depth: " + _maxDepth + ", nodes count: " + _nodesCount);

            CompressedSDT compressedSDT = new CompressedSDT(compress());

            _model._output._treeKey = compressedSDT._key;
            DKV.put(compressedSDT);
            _job.update(1);
            _model.update(_job);
            System.out.println("Tree:");
            System.out.println(Arrays.deepToString(_tree));
        }
        
        private void buildSDTIteratively() {
            _tree = new double[(int) Math.pow(2, _maxDepth + 1)][2];
            Queue<DataFeaturesLimits> limitsQueue = new LinkedList<>();
            limitsQueue.add(getInitialFeaturesLimits());
            for(int nodeIndex = 0; nodeIndex < _tree.length; nodeIndex ++) {
                buildNextNode(limitsQueue, nodeIndex);
            }
        }
        
        

    }


//
//    public void trainSDT() {
//        root = buildSubtree(trainData, getStartFeaturesLimits());
////        System.out.println(root);
//    }

    @Override
    protected Driver trainModelImpl() {
        return new SDTDriver();
    }


    @Override
    public ModelCategory[] can_build() {
        return new ModelCategory[]{
                ModelCategory.Binomial,
//                ModelCategory.Multinomial,
//                                            ModelCategory.Ordinal,
//                ModelCategory.Regression
        };
    }

    @Override
    public boolean isSupervised() {
        return true;
    }

    private static class DataSplit {
        Frame leftSplit;
        Frame rightSplit;
        int featureIndex;
        double threshold;
    }

    /**
     * Not used in new implementation - splits frame
     * @param feature
     * @param threshold
     * @return
     */
    public DataSplit splitData(final int feature, final double threshold) {
        // Define task
        SplitFrameMRTask taskLeftSplit = new SplitFrameMRTask(feature, threshold, 0);
        SplitFrameMRTask taskRightSplit = new SplitFrameMRTask(feature, threshold, 1);
        System.out.println("trainData.types().length: " + _train.types().length);

        DataSplit split = new DataSplit();
        split.featureIndex = feature;
        split.threshold = threshold;
        // Run task
        taskLeftSplit.doAll(_train.types(), _train);
        taskRightSplit.doAll(_train.types(), _train);
        split.leftSplit = taskLeftSplit.outputFrame(
                null, // The output Frame will be stored in DKV and you can access it with this Key, can be null, in case you don't wanted in DKV
                _train.names(),
                _train.domains() // Categorical columns need domain, pass null for Numerical and String columns
        );
        split.rightSplit = taskRightSplit.outputFrame(
                null, // The output Frame will be stored in DKV and you can access it with this Key, can be null, in case you don't wanted in DKV
                _train.names(),
                _train.domains() // Categorical columns need domain, pass null for Numerical and String columns
        );
        return split;
    }

    private Pair<Integer, Integer> countClasses(final DataFeaturesLimits featuresLimits) {
        GetClassCountsMRTask task = new GetClassCountsMRTask(featuresLimits == null
                ? Stream.generate(() -> new double[]{Double.MIN_VALUE, Double.MAX_VALUE}).limit(_train.numCols()).toArray(double[][]::new)
                : featuresLimits.toDoubles());
        task.doAll(_train);

        return new Pair<>(task._count0, task._count1);
    }


    /**
     * Not used now - calculates ratio in the given frame
     * @param data
     * @return
     */
    private Double getZeroRatio(final Frame data) {
        GetClassCountsMRTask task = new GetClassCountsMRTask(null);
        task.doAll(data);

        return task._count0 * 1.0 / (task._count0 + task._count1);
    }

    /**
     * Not used - ratio is now calculated from classes count
     * @param featuresLimits
     * @return
     */
    private Double getZeroRatio(final DataFeaturesLimits featuresLimits) {
        GetClassCountsMRTask task = new GetClassCountsMRTask(featuresLimits == null
                ? Stream.generate(() -> new double[]{Double.MIN_VALUE, Double.MAX_VALUE}).limit(_train.numCols()).toArray(double[][]::new)
                : featuresLimits.toDoubles());
        task.doAll(_train);
//        System.out.println(task._count0 * 1.0 / (task._count0 + task._count1) + " task._count0: " + task._count0 + ", task._count1: " + task._count1);
        return task._count0 * 1.0 / (task._count0 + task._count1);
    }

    public Node getRoot() {
        return _root;
    }


}


// todo - analyze input data to get candidates for threshold
// see random forest with 1 tree
// bin data (go not by step but by count of data)
// data must be solved by histogram because of the size


// todo:
// train random forest on the same data to see the difference
// train sklearn decision tree on the same data
// debug dtree and see how the clustering is used
// more general - n-class classification, regression, gini, different 
// maybe evaluation by depth and hyperparams on validation set inside of training. Read about it