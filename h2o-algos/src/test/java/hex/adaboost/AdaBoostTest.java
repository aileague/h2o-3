package hex.adaboost;

import hex.glm.GLM;
import hex.glm.GLMModel;
import org.apache.commons.io.FileUtils;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.contrib.java.lang.system.EnvironmentVariables;
import org.junit.runner.RunWith;
import water.Scope;
import water.TestUtil;
import water.fvec.Frame;
import water.fvec.TestFrameBuilder;
import water.fvec.Vec;
import water.runner.CloudSize;
import water.runner.H2ORunner;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

@CloudSize(1)
@RunWith(H2ORunner.class)
public class AdaBoostTest extends TestUtil {
    
    public boolean print = true;

    @Rule
    public EnvironmentVariables environmentVariables = new EnvironmentVariables();

    @Before
    public void beforeClass() {
        final File h2oHomeDir = new File(System.getProperty("user.dir")).getParentFile();
        environmentVariables.set("H2O_FILES_SEARCH_PATH", h2oHomeDir.getAbsolutePath());
    }
    
    @Test
    public void testBasicTrain() {
        try {
            Scope.enter();
            Frame train = parseTestFile("smalldata/prostate/prostate.csv");
            String response = "CAPSULE";
            train.toCategoricalCol(response);
            Scope.track(train);
            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._n_estimators = 50;
            p._response_column = response;

            AdaBoost adaBoost = new AdaBoost(p);
            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);
        } finally {
            Scope.exit();
        }
    }

    @Test
    public void testBasicTrainGLM() {
        try {
            Scope.enter();
            Frame train = parseTestFile("smalldata/prostate/prostate.csv");
            Scope.track(train);
            String response = "CAPSULE";
            train.toCategoricalCol(response);
            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._n_estimators = 50;
            p._weak_learner = AdaBoostModel.Algorithm.GLM;
            p._response_column = response;

            AdaBoost adaBoost = new AdaBoost(p);
            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);
        } finally {
            Scope.exit();
        }
    }    

    @Test
    public void testBasicTrainGLMWeakLerner() {
        try {
            Scope.enter();
            Frame train = Scope.track(parseTestFile("smalldata/prostate/prostate.csv"));
            String response = "CAPSULE";
            train.toCategoricalCol(response);
            GLMModel.GLMParameters p = new GLMModel.GLMParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._response_column = response;

            GLM adaBoost = new GLM(p);
            GLMModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);
            Frame score = adaBoostModel.score(train);
            Scope.track(score);
        } finally {
            Scope.exit();
        }
    }    

    @Test
    public void testBasicTrainLarge() {
        try {
            Scope.enter();
            Frame train = parseTestFile("bigdata/laptop/creditcardfraud/creditcardfraud.csv");
            Scope.track(train);
            String response = "Class";
            train.toCategoricalCol(response);
            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._n_estimators = 50;
            p._response_column = response;

            AdaBoost adaBoost = new AdaBoost(p);
            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);
        } finally {
            Scope.exit();
        }
    }

    @Test
    public void testBasicTrainAndScore() {
        try {
            Scope.enter();
            Frame train = parseTestFile("smalldata/prostate/prostate.csv");
            Scope.track(train);
            String response = "CAPSULE";
            train.toCategoricalCol(response);
            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._n_estimators = 50;
            p._response_column = response;

            AdaBoost adaBoost = new AdaBoost(p);
            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);

            System.out.println("train.toTwoDimTable() = " + train.toTwoDimTable());
            
            Frame score = adaBoostModel.score(train);
            Scope.track(score);
            toCSV(score, "../prostatescore.csv");
//            Frame scoreOriginal = Scope.track(parseTestFile("../prostatescore_original.csv"));
//            assertFrameEquals(scoreOriginal, score, 0);
        } finally {
            Scope.exit();
        }
    }

//    @Test
//    public void testBasicTrainAndScoreGLM() {
//        try {
//            Scope.enter();
//            Frame train = Scope.track(parseTestFile("smalldata/prostate/prostate.csv"));
//            Frame test = Scope.track(parseTestFile("smalldata/prostate/prostate.csv"));
//            String response = "CAPSULE";
//            train.toCategoricalCol(response);
//            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
//            p._train = train._key;
//            p._seed = 0xDECAF;
//            p._n_estimators = 2;
//            p._weak_learner = AdaBoostModel.Algorithm.GLM;
//            p._response_column = response;
//
//            AdaBoost adaBoost = new AdaBoost(p);
//            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
//            Scope.track_generic(adaBoostModel);
//            assertNotNull(adaBoostModel);
//
//            Frame score = adaBoostModel.score(test);
//            Scope.track(score);
//            toCSV(score, "../prostatescoreglm.csv");
//        } finally {
//            Scope.exit();
//        }
//    }

    @Test
    public void testBasicTrainAndScoreLarge() {
        try {
            Scope.enter();
            Frame train = parseTestFile("bigdata/laptop/creditcardfraud/creditcardfraud.csv");
            Scope.track(train);
            String response = "Class";
            train.toCategoricalCol(response);
            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._n_estimators = 50;
            p._response_column = response;

            AdaBoost adaBoost = new AdaBoost(p);
            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);

            Frame score = adaBoostModel.score(train);
            Scope.track(score);
            toCSV(score, "../creditcardfraudscore.csv");
        } finally {
            Scope.exit();
        }
    }    

    @Test
    public void testBasicTrainAirlines() {
        try {
            Scope.enter();
            Frame train = parseTestFile("smalldata/testng/airlines_train_preprocessed.csv");
            Scope.track(train);
            Frame test = parseTestFile("smalldata/testng/airlines_test_preprocessed.csv");
            Scope.track(test);
            String response = "IsDepDelayed";
            train.toCategoricalCol(response);
            test.toCategoricalCol(response);
            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._n_estimators = 50;
            p._response_column = response;

            AdaBoost adaBoost = new AdaBoost(p);
            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);

            Frame score = adaBoostModel.score(test);
            Scope.track(score);
            toCSV(score, "../airlinesscore.csv");
        } finally {
            Scope.exit();
        }
    }

    @Test
    public void testBasicTrainHiggs() {
        try {
            Scope.enter();
            Frame train = parseTestFile("smalldata/higgs/higgs_train_5k.csv");
            Scope.track(train);
            Frame test = parseTestFile("smalldata/higgs/higgs_test_5k.csv");
            Scope.track(test);
            String response = "response";
            train.toCategoricalCol(response);
            test.toCategoricalCol(response);
            AdaBoostModel.AdaBoostParameters p = new AdaBoostModel.AdaBoostParameters();
            p._train = train._key;
            p._seed = 0xDECAF;
            p._n_estimators = 50;
            p._response_column = response;

            AdaBoost adaBoost = new AdaBoost(p);
            AdaBoostModel adaBoostModel = adaBoost.trainModel().get();
            Scope.track_generic(adaBoostModel);
            assertNotNull(adaBoostModel);

            Frame score = adaBoostModel.score(test);
            Scope.track(score);
            toCSV(score, "../higgsscore.csv");
        } finally {
            Scope.exit();
        }
    }
    
    private void toCSV(Frame frame, String filename) {
        if (print) {
            File targetFile = new File(filename);
            try {
                FileUtils.copyInputStreamToFile(frame.toCSV(new Frame.CSVStreamParams()), targetFile);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    
    @Test
    public void testCountWe() {
        Scope.enter();
        try {
            Frame train = new TestFrameBuilder()
                    .withVecTypes(Vec.T_NUM, Vec.T_CAT, Vec.T_CAT)
                    .withDataForCol(0, ard(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
                    .withDataForCol(1, ar("0", "0", "0", "0", "0", "1", "1", "1", "1", "1"))
                    .withDataForCol(2, ar("1", "1", "1", "1", "1", "0", "0", "0", "0", "0"))
                    .build();
            train = ensureDistributed(train);
            Scope.track(train);

            CountWeTask countWeTask = new CountWeTask().doAll(train);
            assertEquals("Sum of weights is not correct",10, countWeTask.W, 0);
            assertEquals("Sum of error weights is not correct",10, countWeTask.We, 0);
        } finally {
            Scope.exit();
        }
    }

    @Test
    public void testUpdateWeights() {
        Scope.enter();
        try {
            Frame train = new TestFrameBuilder()
                    .withVecTypes(Vec.T_NUM, Vec.T_CAT, Vec.T_CAT)
                    .withDataForCol(0, ard(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
                    .withDataForCol(1, ar("1", "0", "0", "0", "0", "1", "1", "1", "1", "1"))
                    .withDataForCol(2, ar("1", "1", "1", "1", "1", "0", "0", "0", "0", "0"))
                    .build();
            train = ensureDistributed(train);
            Scope.track(train);

            double alpha = 2;
            UpdateWeightsTask updateWeightsTask = new UpdateWeightsTask(alpha);
            updateWeightsTask.doAll(train);

            Vec weightsExpected =  Vec.makeCon(Math.exp(alpha),train.numRows());
            weightsExpected.set(0, Math.exp(-alpha));
            System.out.println("weights = ");
            System.out.println(new Frame(train.vec(0)).toTwoDimTable(0, (int) train.numRows(), false));
            assertVecEquals("Weights are not correctly updated", weightsExpected, train.vec(0),0);
        } finally {
            Scope.exit();
        }
    }    
}
