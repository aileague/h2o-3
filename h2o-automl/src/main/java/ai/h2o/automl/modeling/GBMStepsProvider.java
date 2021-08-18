package ai.h2o.automl.modeling;

import ai.h2o.automl.*;
import ai.h2o.automl.ModelSelectionStrategies.KeepBestN;
import ai.h2o.automl.events.EventLogEntry;
import hex.Model;
import hex.grid.Grid;
import hex.grid.HyperSpaceSearchCriteria;
import hex.tree.SharedTreeModel;
import hex.tree.gbm.GBMModel;
import hex.tree.gbm.GBMModel.GBMParameters;
import water.Job;
import water.Key;

import java.util.*;

import static ai.h2o.automl.ModelingStep.GridStep.DEFAULT_GRID_TRAINING_WEIGHT;
import static ai.h2o.automl.ModelingStep.ModelStep.DEFAULT_MODEL_TRAINING_WEIGHT;

public class GBMStepsProvider
        implements ModelingStepsProvider<GBMStepsProvider.GBMSteps>
                 , ModelParametersProvider<GBMParameters> {

    public static class GBMSteps extends ModelingSteps {

        static GBMParameters prepareModelParameters() {
            GBMParameters params = new GBMParameters();
            params._score_tree_interval = 5;
            params._histogram_type = SharedTreeModel.SharedTreeParameters.HistogramType.AUTO;
            return params;
        }

        static abstract class GBMModelStep extends ModelingStep.ModelStep<GBMModel> {

            GBMModelStep(String id, int weight, int priorityGroup, AutoML autoML) {
                super(Algo.GBM, id, weight, priorityGroup, autoML);
            }

            protected GBMParameters prepareModelParameters() {
                GBMParameters params = GBMSteps.prepareModelParameters();
                params._ntrees = 10000;
                params._sample_rate = 0.8;
                params._col_sample_rate = 0.8;
                params._col_sample_rate_per_tree = 0.8;
                return params;
            }
        }

        static abstract class GBMGridStep extends ModelingStep.GridStep<GBMModel> {
            public GBMGridStep(String id, int weight, int priorityGroup, AutoML autoML) {
                super(Algo.GBM, id, weight, priorityGroup,autoML);
            }

            protected GBMParameters prepareModelParameters() {
                GBMParameters params = GBMSteps.prepareModelParameters();
                params._ntrees = 10000;
                return params;
            }
        }

        static abstract class GBMExploitationStep extends ModelingStep.SelectionStep<GBMModel> {

            protected GBMModel getBestGBM() {
                for (Model model : getTrainedModels()) {
                    if (model instanceof GBMModel) {
                        return (GBMModel) model;
                    }
                }
                return null;
            }

            @Override
            protected boolean canRun() {
                return super.canRun() && getBestGBM() != null;
            }
            public GBMExploitationStep(String id, int weight, int priorityGroup, AutoML autoML) {
                super(Algo.GBM, id, weight, priorityGroup, autoML);
            }
        }



        private final ModelingStep[] defaults = new GBMModelStep[] {
                new GBMModelStep("def_1", DEFAULT_MODEL_TRAINING_WEIGHT, 3, aml()) {
                    @Override
                    protected GBMParameters prepareModelParameters() {
                        GBMParameters params = super.prepareModelParameters();
                        params._max_depth = 6;
                        params._min_rows = 1;
                        return params;
                    }
                },
                new GBMModelStep("def_2", DEFAULT_MODEL_TRAINING_WEIGHT, 2, aml()) {
                    @Override
                    protected GBMParameters prepareModelParameters() {
                        GBMParameters params = super.prepareModelParameters();
                        params._max_depth = 7;
                        params._min_rows = 10;
                        return params;
                    }
                },
                new GBMModelStep("def_3", DEFAULT_MODEL_TRAINING_WEIGHT, 2, aml()) {
                    @Override
                    protected GBMParameters prepareModelParameters() {
                        GBMParameters params = super.prepareModelParameters();
                        params._max_depth = 8;
                        params._min_rows = 10;
                        return params;
                    }
                },
                new GBMModelStep("def_4", DEFAULT_MODEL_TRAINING_WEIGHT, 2, aml()) {
                    @Override
                    protected GBMParameters prepareModelParameters() {
                        GBMParameters params = super.prepareModelParameters();
                        params._max_depth = 10;
                        params._min_rows = 10;
                        return params;
                    }
                },
                new GBMModelStep("def_5", DEFAULT_MODEL_TRAINING_WEIGHT, 1, aml()) {
                    @Override
                    protected GBMParameters prepareModelParameters() {
                        GBMParameters params = super.prepareModelParameters();
                        params._max_depth = 15;
                        params._min_rows = 100;
                        return params;
                    }
                },
        };

        static class DefaultGBMGridStep extends GBMGridStep {
            public DefaultGBMGridStep(String id, int weight, int priorityGroup, AutoML autoML) {
                super(id, weight, priorityGroup, autoML);
            }

            @Override
            protected Map<String, Object[]> prepareSearchParameters() {
                Map<String, Object[]> searchParams = new HashMap<>();
                searchParams.put("_max_depth", new Integer[]{3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17});
                searchParams.put("_min_rows", new Integer[]{1, 5, 10, 15, 30, 100});
//                        searchParams.put("_learn_rate", new Double[]{0.001, 0.005, 0.008, 0.01, 0.05, 0.08, 0.1, 0.5, 0.8});
                searchParams.put("_sample_rate", new Double[]{0.50, 0.60, 0.70, 0.80, 0.90, 1.00});
                searchParams.put("_col_sample_rate", new Double[]{ 0.4, 0.7, 1.0});
                searchParams.put("_col_sample_rate_per_tree", new Double[]{ 0.4, 0.7, 1.0});
                searchParams.put("_min_split_improvement", new Double[]{1e-4, 1e-5});
                return searchParams;
            }
        }
        
        private final ModelingStep[] grids = new GBMGridStep[] {
                new DefaultGBMGridStep("grid_1", 2*DEFAULT_GRID_TRAINING_WEIGHT, 4, aml()) {
                    @Override
                    protected Job<Grid> startJob() {
                        Job<Grid> job = super.startJob();
                        getResumableResultKeys().put(_algo+"_grid_1", job._result);
                        return job;
                    }
                },
                new DefaultGBMGridStep("grid_1_end", DEFAULT_GRID_TRAINING_WEIGHT, 100, aml()) {
                    @Override
                    protected void setSearchCriteria(HyperSpaceSearchCriteria.RandomDiscreteValueSearchCriteria searchCriteria, Model.Parameters baseParms) {
                        super.setSearchCriteria(searchCriteria, baseParms);
                        searchCriteria.set_stopping_rounds(0);
                    }

                    @Override
                    @SuppressWarnings("unchecked")
                    protected Job<Grid> startJob() {
                        Key<Grid> resumedGrid = getResumableResultKeys().get(_algo+"_grid_1");
                        return hyperparameterSearch(resumedGrid, prepareModelParameters(), prepareSearchParameters());
                    }
                }
        };


        private final ModelingStep[] exploitation = new ModelingStep[] {
                new GBMExploitationStep("lr_annealing", DEFAULT_MODEL_TRAINING_WEIGHT, 6, aml()) {

                    Key<Models> resultKey = null;

                    @Override
                    protected Job<Models> startTraining(Key result, double maxRuntimeSecs) {
                        resultKey = result;
                        GBMModel bestGBM = getBestGBM();
                        aml().eventLog().info(EventLogEntry.Stage.ModelSelection, "Retraining best GBM with learning rate annealing: "+bestGBM._key);
                        GBMParameters params = (GBMParameters) bestGBM._parms.clone();
                        params._ntrees = 10000; // reset ntrees (we'll need more for this fine tuning)
                        params._max_runtime_secs = 0; // reset max runtime
                        params._learn_rate_annealing = 0.99;
                        initTimeConstraints(params, maxRuntimeSecs);
                        setStoppingCriteria(params, new GBMParameters());
                        return asModelsJob(startModel(Key.make(result+"_model"), params), result);
                    }

                    @Override
                    protected ModelSelectionStrategy getSelectionStrategy() {
                        return (originalModels, newModels) ->
                                new KeepBestN<>(1, () -> makeTmpLeaderboard(Objects.toString(resultKey, _algo+"_"+_id)))
                                        .select(new Key[] { getBestGBM()._key }, newModels);
                    }
                }
        };

        public GBMSteps(AutoML autoML) {
            super(autoML);
        }

        @Override
        protected ModelingStep[] getDefaultModels() {
            return defaults;
        }

        @Override
        protected ModelingStep[] getGrids() {
            return grids;
        }

        @Override
        protected ModelingStep[] getExploitation() {
            return exploitation;
        }
    }

    @Override
    public String getName() {
        return Algo.GBM.name();
    }

    @Override
    public GBMSteps newInstance(AutoML aml) {
        return new GBMSteps(aml);
    }

    @Override
    public GBMParameters newDefaultParameters() {
        return new GBMParameters();
    }
}

