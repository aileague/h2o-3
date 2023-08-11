package hex.glm;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import water.DKV;
import water.Scope;
import water.TestUtil;
import water.fvec.Frame;
import water.fvec.TestFrameBuilder;
import water.runner.CloudSize;
import water.runner.H2ORunner;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static hex.glm.GLMModel.GLMParameters.Family.gaussian;
import static hex.glm.GLMModel.GLMParameters.Solver.IRLSM;
import static water.fvec.Vec.T_NUM;
import static water.fvec.Vec.T_STR;

@RunWith(H2ORunner.class)
@CloudSize(1)
public class GLMConstrainedTest extends TestUtil {
  public static final double EPS = 1e-6;
  Frame _betaConstraint1;
  Frame _linearConstraint1;
  Frame _linearConstraint2;
  List<String> _coeffNames1;
  String[][] _betaConstraintNames1;
  double[][] _betaConstraintVal1;
  String[][] _equalityNames1;
  double[][] _equalityValues1;
  String[][] _lessThanNames1;
  double[][] _lessThanValues1;
  
  @Before
  public void createConstraintFrames() {
    Scope.enter();
    Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
    int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
    for (int colInd : catCol)
      train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
    DKV.put(train);
    GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
    params._standardize=true;
    params._response_column = "C21";
    params._max_iterations = 0;
    params._train = train._key;
    GLMModel glm = new GLM(params).trainModel().get();
    _coeffNames1 = Stream.of(glm._output._coefficient_names).collect(Collectors.toList()); ;
    int coefLen = _coeffNames1.size()-1;
    // there are 5 constraints in the beta constraints: 1.0 <= beta0 <= 10.0, -1.0 <= beta1, betacoefLen <= 8.0,
    // 0 <= betacoefLen-1 <= 2.0, 0.1 == betacoefLen-2 == 0.1.  This will be translated into the following standard
    // form: 1.0-beta0 <= 0; beta0-10.0 <= 0, -1.0 -beta1 <= 0, betacoefLen-8.0 <= 0, 0-betacoefLen-1 <= 0,
    // betacoefLen-1-2.0 <= 0, betacoefLen-2-0.1 == 0.
    _betaConstraint1 =
            new TestFrameBuilder()
                    .withColNames("names", "lower_bounds", "upper_bounds")
                    .withVecTypes(T_STR, T_NUM, T_NUM)
                    .withDataForCol(0, new String[] {_coeffNames1.get(0), _coeffNames1.get(1),
                            _coeffNames1.get(coefLen-3), _coeffNames1.get(coefLen-2), _coeffNames1.get(coefLen-1)})
                    .withDataForCol(1, new double [] {1.0, -1.0, Double.NEGATIVE_INFINITY, 0.0, 0.1})
                    .withDataForCol(2, new double[] {10.0, Double.POSITIVE_INFINITY, 8.0, 2.0, 0.1}).build();
    Scope.track(_betaConstraint1);
    // there are two constraints in the linear_constraints, the first one is 2*beta_0+0.5*beta_5 -1<= 0, the second 
    // one is 0.5*beta_36-1.5*beta_38 == 0
    _linearConstraint1 = new TestFrameBuilder()
            .withColNames("names", "values", "types", "constraint_numbers")
            .withVecTypes(T_STR, T_NUM, T_STR, T_NUM)
            .withDataForCol(0, new String[] {_coeffNames1.get(0), _coeffNames1.get(5), "constant",
                    _coeffNames1.get(36), _coeffNames1.get(38)})
            .withDataForCol(1, new double [] {2, 0.5, -1, 0.5, -1.5})
            .withDataForCol(2, new String[] {"lessthanequal", "lessthanequal", "lessthanequal", "equal", "equal"})
            .withDataForCol(3, new int[]{0,0,0,1,1}).build();
    Scope.track(_linearConstraint1);
    // form correct constraints names and values:
    _betaConstraintNames1 = new String[][]{{_coeffNames1.get(0), "constant"}, {_coeffNames1.get(0), "constant"},
            {_coeffNames1.get(1), "constant"}, {_coeffNames1.get(coefLen-3), "constant"},
            {_coeffNames1.get(coefLen-2), "constant"}, {_coeffNames1.get(coefLen-2), "constant"},
            {_coeffNames1.get(coefLen-1), "constant"}};
    _betaConstraintVal1 = new double[][]{{-1,1}, {1,-10}, {-1,-1},
            {1.0,-8*train.vec(_coeffNames1.get(coefLen-3)).sigma()},
            {-1,0*train.vec(_coeffNames1.get(coefLen-2)).sigma()},
            {1,-2.0*train.vec(_coeffNames1.get(coefLen-2)).sigma()}, {1,-0.1*train.vec(_coeffNames1.get(coefLen-1)).sigma()}};

    _equalityNames1 = new String[][]{{_coeffNames1.get(36), _coeffNames1.get(38), "constant"}};
    _equalityValues1 = new double[][]{{0.5/train.vec(_coeffNames1.get(36)).sigma(),
            -1.5/train.vec(_coeffNames1.get(38)).sigma(), 0.0}};

    _lessThanNames1 = new String[][]{{_coeffNames1.get(0), _coeffNames1.get(5), "constant"}};
    _lessThanValues1 = new double[][]{{2, 0.5, -1}};

    // Constraints in the linear_constraints, 
    // a. -0.3*beta_0+0.5*beta_1+1*beta_3-3 <= 0; 
    // b. 3*beta_8-4*beta_36+0.5*beta_37 <= 0, 
    // c.0.1*beta_38-0.2*beta_39==0, 
    // d. 2*beta_40-0.1*beta_41-0.4*beta_42+0.8 <= 0;
    // e. 0.1*beta_4-0.5*beta_5+0.7*beta_6-1.1 == 0; 
    // f. 2*beta_6+0.5*beta_43-0.3*beta_7 == 0 
    // g. 0.5*beta_36-1.5*beta_38-0.3 == 0    
    _linearConstraint2 = new TestFrameBuilder()
            .withColNames("names", "values", "types", "constraint_numbers")
            .withVecTypes(T_STR, T_NUM, T_STR, T_NUM)
            .withDataForCol(0, new String[] {_coeffNames1.get(0), _coeffNames1.get(1), _coeffNames1.get(3),
                    "constant", _coeffNames1.get(8), _coeffNames1.get(36), _coeffNames1.get(37), _coeffNames1.get(38),
                    _coeffNames1.get(39), _coeffNames1.get(40), _coeffNames1.get(41), _coeffNames1.get(42), "constant",
                    _coeffNames1.get(4), _coeffNames1.get(5), _coeffNames1.get(6), "constant", _coeffNames1.get(6),
                    _coeffNames1.get(43), _coeffNames1.get(7), _coeffNames1.get(36), _coeffNames1.get(38), "constant"})
            .withDataForCol(1, new double [] {-0.3, 0.5, 1.0, -3.0, 3, -4, 0.5, 0.1, -0.2, 2.0, -0.1, -0.4,
                    0.8, 0.1, -0.5, 0.7, -1.1, 2.0, 0.5, -0.3, 0.5, -1.5, -0.3})
            .withDataForCol(2, new String[] {"lessthanequal", "lessthanequal", "lessthanequal", "lessthanequal",
                    "lessthanequal", "lessthanequal", "lessthanequal", "equal", "equal", "lessthanequal",
                    "lessthanequal", "lessthanequal", "lessthanequal", "equal", "equal", "equal", "equal", "equal",
                    "equal", "equal", "equal", "equal", "equal"})
            .withDataForCol(3, new int[]{0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6,
                    6}).build();
    Scope.track(_linearConstraint2);
    Scope.untrack(train);

    glm.delete();
  }
  
  @After
  public void teardown() {
    Scope.exit();
  }
  
  // make sure correct constraint matrix is formed after extracting constraints from beta constraints and linear
  // constraints
  @Test
  public void testBetaLinearConstraintMatrix() {
    Scope.enter();
    try {
      Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
      int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
      for (int colInd : catCol)
        train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
      DKV.put(train);
      GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
      params._standardize = true;
      params._response_column = "C21";
      params._solver = IRLSM;
      params._train = train._key;
      // build the beta_constraints
      Frame beta_constraints = _betaConstraint1;
      Frame linear_constraints = _linearConstraint1;
      params._max_iterations = 1;
      params._expose_constraints = true;
      params._beta_constraints = beta_constraints._key;
      params._linear_constraints = linear_constraints._key;
      GLMModel glm2 = new GLM(params).trainModel().get();
      Scope.track_generic(glm2);
      // check that constraint matrix is extracted correctly
      List<String> constraintNames = Arrays.stream(glm2._output._constraintCoefficientNames).collect(Collectors.toList());
      double[][] initConstraintMatrix = glm2._output._initConstraintMatrix;
      // check rows from beta constraints
      int rowIndex = 0;
      assertCorrectConstraintMatrix(constraintNames, initConstraintMatrix, _betaConstraintNames1, _betaConstraintVal1, rowIndex);
      // check rows from linear constraints with equality
      rowIndex += _betaConstraintNames1.length;
      assertCorrectConstraintMatrix(constraintNames, initConstraintMatrix, _equalityNames1, _equalityValues1, rowIndex);
      // check row from linear contraints with lessThanEqualTo
      rowIndex += _equalityNames1.length;
      assertCorrectConstraintMatrix(constraintNames, initConstraintMatrix, _lessThanNames1, _lessThanValues1, rowIndex);
    } finally {
      Scope.exit();
    }
  }
  
  public void assertCorrectConstraintMatrix(List<String> constraintNames, double[][] constraintMatrix,
                                            String[][] origNames, double[][] originalValues, int rowStart) {
    int numConstraints = origNames.length;
    for (int index=0; index<numConstraints; index++) {
      int rowIndex = index+rowStart;
      String[] constNames = origNames[index];
      double[] constValues = originalValues[index];
      int numNames = constNames.length;
      for (int index2=0; index2<numNames; index2++) {
        int cNamesIndex = constraintNames.indexOf(constNames[index2]);
        assert Math.abs(constraintMatrix[rowIndex][cNamesIndex]-constValues[index2])< EPS : 
                "Expected valud: "+constValues[index2]+" for constraint "+constNames[index2]+" but actual: "
                        +constraintMatrix[rowIndex][cNamesIndex];
      }
    }
  }
  
  // make sure we can get coefficient names without building a GLM model.  We compare the coefficient names
  // obtained without building a model and with building a model.  They should be the same.
  @Test
  public void testCoefficientNames() {
    Scope.enter();
    try {
      Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
      int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
      for (int colInd : catCol)
        train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
      DKV.put(train);
      GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
      params._standardize=true;
      params._response_column = "C21";
      params._max_iterations = 0;
      params._train = train._key;
      GLMModel glm = new GLM(params).trainModel().get();
      Scope.track_generic(glm);
      List<String> coeffNames = Stream.of(glm._output._coefficient_names).collect(Collectors.toList()); ;
  
      params._max_iterations = 1;
      GLMModel glm2 = new GLM(params).trainModel().get();
      Scope.track_generic(glm2);
      String[] coeffNames2 = glm2.coefficients().keySet().toArray(new String[0]);
      
      for (String oneName : coeffNames2)
        assert coeffNames.contains(oneName);
    } finally {
      Scope.exit();
    }
  }

  // test constraints specified in beta_constraint and linear constraints and extracted correctly with 
  // standardization.
  @Test
  public void testConstraintsInBetaLinearStandard() {
    Scope.enter();
    try {
      Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
      int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
      for (int colInd : catCol)
        train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
      DKV.put(train);
      GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
      params._standardize = true;
      params._response_column = "C21";
      params._solver = IRLSM;
      params._train = train._key;
      List<String> coeffNames = _coeffNames1;
      // build the beta_constraints
      int coefLen = coeffNames.size()-1;
      Frame beta_constraints = _betaConstraint1;
      Frame linear_constraints = _linearConstraint1;
      params._max_iterations = 1;
      params._expose_constraints = true;
      params._beta_constraints = beta_constraints._key;
      params._linear_constraints = linear_constraints._key;
      GLMModel glm2 = new GLM(params).trainModel().get();
      Scope.track_generic(glm2);
      // check constraints from betaConstraints are extracted properly
      ConstrainedGLMUtils.LinearConstraints[] fromBetaConstraint = glm2._output._fromBetaConstraints;
      assert _betaConstraintNames1.length == fromBetaConstraint.length : "Expected constraint length: "+_betaConstraintVal1.length+" but" +
              " actual is "+fromBetaConstraint.length;
      assertCorrectConstraintContent(_betaConstraintNames1, _betaConstraintVal1, fromBetaConstraint);
      // check constraints from linear constraints are extracted properly
      // check equality constraint
      assertCorrectConstraintContent(_equalityNames1, _equalityValues1, glm2._output._equalityConstraints);
      // check lessThanEqual to constraint
      assertCorrectConstraintContent(_lessThanNames1, _lessThanValues1, glm2._output._lessThanEqualToConstraints);
    } finally {
      Scope.exit();
    }
  }


  // test constraints specified in beta_constraint and linear constraints and extracted correctly without 
  // standardization.
  @Test
  public void testConstraintsInBetaLinear() {
    Scope.enter();
    try {
      Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
      int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
      for (int colInd : catCol)
        train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
      DKV.put(train);
      List<String> coeffNames = _coeffNames1;
      // build the beta_constraints
      int coefLen = coeffNames.size()-1;
      Frame beta_constraints = _betaConstraint1;
      // there are two constraints in the linear_constraints, the first one is 2*beta_0+0.5*beta_5 -1<= 0, the second 
      // one is 0.5*beta_36-1.5*beta_38 == 0
      Frame linear_constraints = _linearConstraint1;
      GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
      params._standardize=false;
      params._response_column = "C21";
      params._max_iterations = 0;
      params._solver = IRLSM;
      params._train = train._key;
      params._max_iterations = 1;
      params._expose_constraints = true;
      params._beta_constraints = beta_constraints._key;
      params._linear_constraints = linear_constraints._key;
      GLMModel glm2 = new GLM(params).trainModel().get();
      Scope.track_generic(glm2);
      // check constraints from betaConstraints are extracted properly
      String[][] betaConstraintNames = new String[][]{{coeffNames.get(0), "constant"}, {coeffNames.get(0), "constant"}, 
              {coeffNames.get(1), "constant"}, {coeffNames.get(coefLen-3), "constant"}, 
              {coeffNames.get(coefLen-2), "constant"}, {coeffNames.get(coefLen-2), "constant"}, 
              {coeffNames.get(coefLen-1), "constant"}};
      double[][] betaConstraintVal = new double[][]{{-1,1}, {1,-10}, {-1,-1}, {1,-8}, {-1,0}, {1,-2}, {1,-0.1}};
      ConstrainedGLMUtils.LinearConstraints[] fromBetaConstraint = glm2._output._fromBetaConstraints;
      assert betaConstraintNames.length == fromBetaConstraint.length : "Expected constraint length: "+betaConstraintVal.length+" but" +
              " actual is "+fromBetaConstraint.length;
      assertCorrectConstraintContent(betaConstraintNames, betaConstraintVal, fromBetaConstraint);
      // check constraints from linear constraints are extracted properly
      // check equality constraint
      String[][] equalityNames = new String[][]{{coeffNames.get(36), coeffNames.get(38), "constant"}};
      double[][] equalityValues = new double[][]{{0.5, -1.5, 0.0}};
      assertCorrectConstraintContent(equalityNames, equalityValues, glm2._output._equalityConstraints);
      // check lessThanEqual to constraint
      String[][] lessThanNames = new String[][]{{coeffNames.get(0), coeffNames.get(5), "constant"}};
      double[][] lessThanValues = new double[][]{{2, 0.5, -1}};
      assertCorrectConstraintContent(lessThanNames, lessThanValues, glm2._output._lessThanEqualToConstraints);
    } finally {
      Scope.exit();
    }
  }
  

  // test constraints specified only in linear_constraint and extracted correctly without standardization
  @Test
  public void testConstraintsInLinear() {
    Scope.enter();
    try {
      Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
      int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
      for (int colInd : catCol)
        train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
      DKV.put(train);
      List<String> coeffNames = _coeffNames1;
      Frame linear_constraints = _linearConstraint2;
      
      GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
      params._standardize=false;
      params._response_column = "C21";
      params._max_iterations = 0;
      params._solver = IRLSM;
      params._train = train._key;
      params._max_iterations = 1;
      params._expose_constraints = true;
      params._linear_constraints = linear_constraints._key;
      GLMModel glm2 = new GLM(params).trainModel().get();
      Scope.track_generic(glm2);
      // check constraints from linear constraints are extracted properly
      // check equality constraint
      String[][] equalityNames = new String[][]{{coeffNames.get(38), coeffNames.get(39), "constant"}, 
              {coeffNames.get(4), coeffNames.get(5), coeffNames.get(6), "constant"},
              {coeffNames.get(6), coeffNames.get(43), coeffNames.get(7), "constant"},
              {coeffNames.get(36), coeffNames.get(38), "constant"}};
      double[][] equalityValues = new double[][]{{0.1, -0.2, 0.0}, {0.1, -0.5, 0.7, -1.1}, {2, 0.5, -0.3, 0.0}, 
              {0.5, -1.5, -0.3}};
      assertCorrectConstraintContent(equalityNames, equalityValues, glm2._output._equalityConstraints);
      // check lessThanEqual to constraint
      String[][] lessThanNames = new String[][]{{coeffNames.get(0), coeffNames.get(1), coeffNames.get(3),
              "constant"}, {coeffNames.get(8), coeffNames.get(36), coeffNames.get(37), "constant"}, 
              {coeffNames.get(40), coeffNames.get(41), coeffNames.get(42), "constant"}};
      double[][] lessThanValues = new double[][]{{-0.3, 0.5, 1, -3}, {3, -4, 0.5, 0.0}, {2, -0.1, -0.4, 0.8}};
      assertCorrectConstraintContent(lessThanNames, lessThanValues, glm2._output._lessThanEqualToConstraints);
    } finally {
      Scope.exit();
    }
  }

  // test constraints specified only in linear_constraint and extracted correctly with standardization
  @Test
  public void testConstraintsInLinearStandard() {
    Scope.enter();
    try {
      Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
      int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
      for (int colInd : catCol)
        train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
      DKV.put(train);
      GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
      params._standardize = true;
      params._response_column = "C21";
      params._max_iterations = 0;
      params._solver = IRLSM;
      params._train = train._key;
      GLMModel glm = new GLM(params).trainModel().get();
      Scope.track_generic(glm);
      List<String> coeffNames = Stream.of(glm._output._coefficient_names).collect(Collectors.toList()); ;
      // build the beta_constraints
      int coefLen = coeffNames.size()-1;

      Frame linear_constraints = _linearConstraint2;
      params._max_iterations = 1;
      params._expose_constraints = true;
      params._linear_constraints = linear_constraints._key;
      GLMModel glm2 = new GLM(params).trainModel().get();
      Scope.track_generic(glm2);
      // check constraints from linear constraints are extracted properly
      // check equality constraint
      String[][] equalityNames = new String[][]{{coeffNames.get(38), coeffNames.get(39), "constant"},
              {coeffNames.get(4), coeffNames.get(5), coeffNames.get(6), "constant"},
              {coeffNames.get(6), coeffNames.get(43), coeffNames.get(7), "constant"},
              {coeffNames.get(36), coeffNames.get(38), "constant"}};
      double[][] equalityValues = new double[][]{{0.1/train.vec(coeffNames.get(38)).sigma(), 
              -0.2/train.vec(coeffNames.get(39)).sigma(), 0.0}, {0.1, -0.5, 0.7, -1.1}, 
              {2, 0.5/train.vec(coeffNames.get(43)).sigma(), -0.3, 0.0},
              {0.5/train.vec(coeffNames.get(36)).sigma(), -1.5/train.vec(coeffNames.get(38)).sigma(), -0.3}};
      assertCorrectConstraintContent(equalityNames, equalityValues, glm2._output._equalityConstraints);
      // check lessThanEqual to constraint
      String[][] lessThanNames = new String[][]{{coeffNames.get(0), coeffNames.get(1), coeffNames.get(3),
              "constant"}, {coeffNames.get(8), coeffNames.get(36), coeffNames.get(37), "constant"},
              {coeffNames.get(40), coeffNames.get(41), coeffNames.get(42), "constant"}};
      double[][] lessThanValues = new double[][]{{-0.3, 0.5, 1, -3}, 
              {3, -4/train.vec(coeffNames.get(36)).sigma(), 0.5/train.vec(coeffNames.get(37)).sigma(), 0.0}, 
              {2/train.vec(coeffNames.get(40)).sigma(), -0.1/train.vec(coeffNames.get(41)).sigma(), 
                      -0.4/train.vec(coeffNames.get(42)).sigma(), 0.8}};
      assertCorrectConstraintContent(lessThanNames, lessThanValues, glm2._output._lessThanEqualToConstraints);
    } finally {
      Scope.exit();
    }
  }

  // Test to make sure the matrix formed from constraints are correct
  @Test
  public void testConstraintsMatrix() {
    Scope.enter();
    try {
      Frame train = parseAndTrackTestFile("smalldata/glm_test/gaussian_20cols_10000Rows.csv");
      int[] catCol = new int[]{0,1,2,3,4,5,6,7,8,9};
      for (int colInd : catCol)
        train.replace((colInd), train.vec(colInd).toCategoricalVec()).remove();
      DKV.put(train);
      GLMModel.GLMParameters params = new GLMModel.GLMParameters(gaussian);
      params._standardize = false;
      params._response_column = "C21";
      params._max_iterations = 0;
      params._solver = IRLSM;
      params._train = train._key;
      GLMModel glm = new GLM(params).trainModel().get();
      Scope.track_generic(glm);
      List<String> coeffNames = Stream.of(glm._output._coefficient_names).collect(Collectors.toList()); ;
      // build the beta_constraints
      int coefLen = coeffNames.size()-1;
      // there are 5 constraints in the beta constraints: 1.0 <= beta0 <= 10.0, -1.0 <= beta1, betacoefLen <= 8.0,
      // 0 <= betacoefLen-1 <= 2.0, 0.1 == betacoefLen-2 == 0.1.  This will be translated into the following standard
      // form: 1.0-beta0 <= 0; beta0-10.0 <= 0, -1.0 -beta1 <= 0, betacoefLen-8.0 <= 0, 0-betacoefLen-1 <= 0,
      // betacoefLen-1-2.0 <= 0, betacoefLen-2-0.1 == 0. 
      Frame beta_constraints =
              new TestFrameBuilder()
                      .withColNames("names", "lower_bounds", "upper_bounds")
                      .withVecTypes(T_STR, T_NUM, T_NUM)
                      .withDataForCol(0, new String[] {coeffNames.get(0), coeffNames.get(1),
                              coeffNames.get(coefLen-3), coeffNames.get(coefLen-2), coeffNames.get(coefLen-1)})
                      .withDataForCol(1, new double [] {1.0, -1.0, Double.NEGATIVE_INFINITY, 0.0, 0.1})
                      .withDataForCol(2, new double[] {10.0, Double.POSITIVE_INFINITY, 8.0, 2.0, 0.1}).build();
      Scope.track(beta_constraints);
      // there are two constraints in the linear_constraints, the first one is 2*beta_0+0.5*beta_5 -1<= 0, the second 
      // one is 0.5*beta_36-1.5*beta_38 == 0
      Frame linear_constraints = new TestFrameBuilder()
              .withColNames("names", "values", "types", "constraint_numbers")
              .withVecTypes(T_STR, T_NUM, T_STR, T_NUM)
              .withDataForCol(0, new String[] {coeffNames.get(0), coeffNames.get(5), "constant",
                      coeffNames.get(36), coeffNames.get(38)})
              .withDataForCol(1, new double [] {2, 0.5, -1, 0.5, -1.5})
              .withDataForCol(2, new String[] {"lessthanequal", "lessthanequal", "lessthanequal", "equal", "equal"})
              .withDataForCol(3, new int[]{0,0,0,1,1}).build();
      params._max_iterations = 1;
      params._expose_constraints = true;
      params._beta_constraints = beta_constraints._key;
      params._linear_constraints = linear_constraints._key;
      GLMModel glm2 = new GLM(params).trainModel().get();
      Scope.track_generic(glm2);
      

    } finally {
      Scope.exit();
    }
  }

  public void assertCorrectConstraintContent(String[][] coefNames, double[][] value,
                                             ConstrainedGLMUtils.LinearConstraints[] consts) {
    assert coefNames.length == consts.length;
    int constLen = consts.length;
    for (int index=0; index<constLen; index++) {
      ConstrainedGLMUtils.LinearConstraints oneConstraint = consts[index];
      Set<String> coefKeys = oneConstraint._constraints.keySet();
      String[] coefName = coefNames[index];
      int entryLen = coefName.length;
      for (int ind = 0; ind < entryLen; ind++) {
        assert coefKeys.contains(coefName[ind]);
        assert Math.abs(value[index][ind] - oneConstraint._constraints.get(coefName[ind])) < EPS : "Expected: "+value[index][ind]+
                ".  Actual: "+oneConstraint._constraints.get(coefName[ind])+".";
      }
    }
  }
}
