package hex.glm;

import Jama.Matrix;
import hex.DataInfo;
import water.DKV;
import water.Iced;
import water.Key;
import water.Scope;
import water.fvec.Frame;
import water.util.IcedHashMap;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class ConstrainedGLMUtils {
  public static class LinearConstraints extends Iced { // store one linear constraint
    IcedHashMap<String, Double> _constraints; // column names, coefficient of constraints
    double _constraintsVal; // contains evaluated constraint values
    
    public LinearConstraints() {
      _constraints = new IcedHashMap<>();
      _constraintsVal = Double.NaN; // represent constraint not evaluated.
    }
  }

  /***
   *
   * This method will extract the constraints specified in beta constraint and combine it with the linear constraints
   * later.  Note that the linear constraints are only accepted in standard form, meaning we only accept the following
   * constraint forms:  2*beta_1-3*beta_4-3 == 0 or 2*beta_1-3*beta_4-3 <= 0.
   * 
   * The beta constraints on the other hand is specified in several forms:
   * 1): -Infinity <= beta <= Infinity: ignored, no constrain here;
   * 2): -Infinity <= beta <= high_val: transformed to beta - high_val <= 0, add to lessThanEqualTo constraint;
   * 3): low_val <= beta <= Infinity: transformed to low_val - beta <= 0, add to lessThanEqualTo constraint;
   * 4): low_val <= beta <= high_val: transformed to two constraints, low_val-beta <= 0, beta-high_val <= 0, add to lessThanEqualTo constraint;
   * 5): val <= beta <= val: transformed to beta-val == 0, add to equalTo constraint.
   * 
   * The newly extracted constraints will be added to fields in state.
   * 
   */
  public static void extractBetaConstraints(ComputationState state, String[] coefNames) {
    GLM.BetaConstraint betaC = state.activeBC();
    List<LinearConstraints> betaConstraints = new ArrayList<>();
    if (betaC._betaLB != null) {
      int numCons = betaC._betaLB.length-1;
      for (int index=0; index<numCons; index++) {
        if (!Double.isInfinite(betaC._betaUB[index]) && (betaC._betaLB[index] == betaC._betaUB[index])) { // equality constraint
          addBCEqualityConstraint(betaConstraints, betaC, coefNames, index);
        } else if (!Double.isInfinite(betaC._betaUB[index]) && !Double.isInfinite(betaC._betaLB[index]) && 
                (betaC._betaLB[index] < betaC._betaUB[index])) { // low < beta < high, generate two lessThanEqualTo constraints
          addBCGreaterThanConstraint(betaConstraints, betaC, coefNames, index);
          addBCLessThanConstraint(betaConstraints, betaC, coefNames, index);
        } else if (Double.isInfinite(betaC._betaUB[index]) && !Double.isInfinite(betaC._betaLB[index])) {  // low < beta < positive infinity
          addBCGreaterThanConstraint(betaConstraints, betaC, coefNames, index);
        } else if (!Double.isInfinite(betaC._betaUB[index]) && Double.isInfinite(betaC._betaLB[index])) { // negative infinity < beta < high
          addBCLessThanConstraint(betaConstraints, betaC, coefNames, index);
        }
      }
    }
    state._fromBetaConstraints = betaConstraints.toArray(new LinearConstraints[0]);
  }

  /***
   * This method will extract the equality constraint and add to equalityC from beta constraint by doing the following
   * transformation: val <= beta <= val: transformed to beta-val == 0, add to equalTo constraint.
   */
  public static void addBCEqualityConstraint(List<LinearConstraints> equalityC, GLM.BetaConstraint betaC,
                                           String[] coefNames, int index) {
    LinearConstraints oneEqualityConstraint = new LinearConstraints();
    oneEqualityConstraint._constraints.put(coefNames[index], 1.0);
    oneEqualityConstraint._constraints.put("constant", -betaC._betaLB[index]);
    equalityC.add(oneEqualityConstraint);
  }

  /***
   * This method will extract the greater than constraint and add to lessThanC from beta constraint by doing the following
   * transformation: low_val <= beta <= Infinity: transformed to low_val - beta <= 0.
   */
  public static void addBCGreaterThanConstraint(List<LinearConstraints> lessThanC, GLM.BetaConstraint betaC,
                                             String[] coefNames, int index) {
    LinearConstraints lessThanEqualToConstraint = new LinearConstraints();
    lessThanEqualToConstraint._constraints.put(coefNames[index], -1.0);
    lessThanEqualToConstraint._constraints.put("constant", betaC._betaLB[index]);
    lessThanC.add(lessThanEqualToConstraint);
  }

  /***
   * This method will extract the less than constraint and add to lessThanC from beta constraint by doing the following
   * transformation: -Infinity <= beta <= high_val: transformed to beta - high_val <= 0.
   */
  public static void addBCLessThanConstraint(List<LinearConstraints> lessThanC, GLM.BetaConstraint betaC,
                                             String[] coefNames, int index) {
    LinearConstraints greaterThanConstraint = new LinearConstraints();
    greaterThanConstraint._constraints.put(coefNames[index], 1.0);
    greaterThanConstraint._constraints.put("constant", -betaC._betaUB[index]);
    lessThanC.add(greaterThanConstraint);  
  }

  /***
   * This method will extract the constraints specified in the Frame with key linearCOnstraintFrameKey.  For example,
   * the following constraints a*beta_1+b*beta_2-c*beta_5 == 0, d*beta_2+e*beta_6-f <= 0 can be specified as the
   * following rows:
   *  names           values            Type            constraint_numbers
   *  beta_1             a              Equal                   0
   *  beta_2             b              Equal                   0
   *  beta_5            -c              Equal                   0
   *  beta_2             d              LessThanEqual           1
   *  beta_6             e              LessThanEqual           1
   *  constant          -f              LessThanEqual           1
   */
  public static void extractLinearConstraints(ComputationState state, Key<Frame> linearConstraintFrameKey, DataInfo dinfo) {
    List<LinearConstraints> equalityC = new ArrayList<>();
    List<LinearConstraints> lessThanEqualToC = new ArrayList<>();
    Frame linearConstraintF = DKV.getGet(linearConstraintFrameKey);
    Scope.track(linearConstraintF);
    List<String> colNamesList = Stream.of(dinfo._adaptedFrame.names()).collect(Collectors.toList());
    List<String> coefNamesList = Stream.of(dinfo.coefNames()).collect(Collectors.toList());
    int numberOfConstraints = linearConstraintF.vec("constraint_numbers").toCategoricalVec().domain().length;
    int numRow = (int) linearConstraintF.numRows();
    List<Integer> rowIndices = IntStream.range(0,numRow).boxed().collect(Collectors.toList());
    String constraintType;
    int rowIndex;
    for (int conInd = 0; conInd < numberOfConstraints; conInd++) {
      if (!rowIndices.isEmpty()) {
        rowIndex = rowIndices.get(0);
        constraintType = linearConstraintF.vec("types").stringAt(rowIndex).toLowerCase();
        if ("equal".equals(constraintType)) {
          extractConstraint(linearConstraintF, rowIndices, equalityC, dinfo, coefNamesList, colNamesList);
        } else if ("lessthanequal".equals(constraintType)) {
          extractConstraint(linearConstraintF, rowIndices, lessThanEqualToC, dinfo, coefNamesList,
                  colNamesList);
        } else {
          throw new IllegalArgumentException("Type of linear constraints can only be Equal to LessThanEqualTo.");
        }
      }
    }
    state.setLinearConstraints(equalityC.toArray(new LinearConstraints[0]), 
            lessThanEqualToC.toArray(new LinearConstraints[0]));
  }
  
  public static void extractConstraint(Frame constraintF, List<Integer> rowIndices, List<LinearConstraints> equalC,
                                       DataInfo dinfo, List<String> coefNames, List<String> colNames) {
    List<Integer> processedRowIndices = new ArrayList<>();
    int constraintNumberFrame = (int) constraintF.vec("constraint_numbers").at(rowIndices.get(0));
    LinearConstraints currentConstraint = new LinearConstraints();
    String constraintType = constraintF.vec("types").stringAt(rowIndices.get(0)).toLowerCase();
    boolean standardize = dinfo._normMul != null;
    boolean constantFound = false;
    for (Integer rowIndex : rowIndices) {
      String coefName = constraintF.vec("names").stringAt(rowIndex);
      String currType = constraintF.vec("types").stringAt(rowIndex).toLowerCase();
      if (!coefNames.contains(coefName) && !"constant".equals(coefName))
        throw new IllegalArgumentException("Coefficient name " + coefName + " is not a valid coefficient name.  It " +
                "be a valid coefficient name or it can be constant");
      if ((int) constraintF.vec("constraint_numbers").at(rowIndex) == constraintNumberFrame) {
        if (!constraintType.equals(currType))
          throw new IllegalArgumentException("Constraint type "+" of the same constraint must be the same but is not." +
                  "  Expected type: "+constraintType+".  Actual type: "+currType);
        if ("constant".equals(coefName))
          constantFound = true;
        processedRowIndices.add(rowIndex);
        // coefNames is valid
        if (standardize && colNames.contains(coefName)) {  // numerical column with standardization
          int colInd = colNames.indexOf(coefName);
          currentConstraint._constraints.put(coefName, constraintF.vec("values").at(rowIndex)*dinfo._normMul[colInd-dinfo._cats]);
        } else {  // categorical column, constant or numerical column without standardization
          currentConstraint._constraints.put(coefName, constraintF.vec("values").at(rowIndex));
        }
      }
    }
    if (!constantFound)
      currentConstraint._constraints.put("constant", 0.0);  // put constant of 0.0
    if (currentConstraint._constraints.size() < 3)
      throw new IllegalArgumentException("linear constraint must have at least two coefficients.  For constraints on" +
              " just one coefficient: "+ constraintF.vec("names").stringAt(0)+", use betaConstraints instead.");
    equalC.add(currentConstraint);
    rowIndices.removeAll(processedRowIndices);
  }
  
  public static double[][] formConstraintMatrix(ComputationState state, List<String> constraintNamesList) {
    // extract coefficient names from constraints
    constraintNamesList.addAll(extractConstraintCoeffs(state));
    // form double matrix
    int numRow = (state._fromBetaConstraints == null ? 0 : state._fromBetaConstraints.length) +
            (state._equalityConstraints == null ? 0 : state._equalityConstraints.length) + 
            (state._lessThanEqualToConstraints == null ? 0 : state._lessThanEqualToConstraints.length);
    double[][] initConstraintMatrix = new double[numRow][constraintNamesList.size()];
    fillConstraintValues(state, constraintNamesList, initConstraintMatrix);
    return initConstraintMatrix;
  }
  
  public static void fillConstraintValues(ComputationState state, List<String> constraintNamesList, double[][] initCMatrix) {
    int rowIndex = 0;
    if (state._fromBetaConstraints != null)
      rowIndex = extractConstraintValues(state._fromBetaConstraints, constraintNamesList, initCMatrix, rowIndex);
    if (state._equalityConstraints != null)
      rowIndex = extractConstraintValues(state._equalityConstraints, constraintNamesList, initCMatrix, rowIndex);
    if (state._lessThanEqualToConstraints != null)
      extractConstraintValues(state._lessThanEqualToConstraints, constraintNamesList, initCMatrix, rowIndex);
  }
  
  public static int extractConstraintValues(LinearConstraints[] constraints, List<String> constraintNamesList, double[][] initCMatrix, int rowIndex) {
    int numConstr = constraints.length;
    for (int index=0; index<numConstr; index++) {
      Set<String> coeffKeys = constraints[index]._constraints.keySet();
      for (String oneKey : coeffKeys) {
        if ( constraintNamesList.contains(oneKey))
          initCMatrix[rowIndex][constraintNamesList.indexOf(oneKey)] = constraints[index]._constraints.get(oneKey);
      }
      rowIndex++;
    }
    return rowIndex;
  }
  
  public static List<String> extractConstraintCoeffs(ComputationState state) {
    List<String> tConstraintCoeffName = new ArrayList<>();
    boolean nonZeroConstant = false;
    if (state._fromBetaConstraints != null)   // extract coefficients constraints
      nonZeroConstant = nonZeroConstant | extractCoeffNames(tConstraintCoeffName, state._fromBetaConstraints);

    if (state._equalityConstraints != null)
      nonZeroConstant = nonZeroConstant | extractCoeffNames(tConstraintCoeffName, state._equalityConstraints);

    if (state._lessThanEqualToConstraints != null)
      nonZeroConstant = nonZeroConstant | extractCoeffNames(tConstraintCoeffName, state._lessThanEqualToConstraints);
    
    // remove duplicates in the constraints names
    Set<String> noDuplicateNames = new HashSet<>(tConstraintCoeffName);
    if (!nonZeroConstant) // no non-Zero constant present
      noDuplicateNames.remove("constant");
    return new ArrayList<>(noDuplicateNames);
  }
  
  public static boolean extractCoeffNames(List<String> coeffList, LinearConstraints[] constraints) {
    int numConst = constraints.length;
    boolean nonZeroConstant = false;
    for (int index=0; index<numConst; index++) {
      Set<String> keys = constraints[index]._constraints.keySet();
      coeffList.addAll(keys);
      if (keys.contains("constant"))
        nonZeroConstant = constraints[index]._constraints.get("constant") != 0.0;
    }
    return nonZeroConstant;
  }
  
  public static List<String> foundRedundantConstraints(ComputationState state, final double[][] initConstraintMatrix) {
    Matrix constMatrix = new Matrix(initConstraintMatrix);
    Matrix constMatrixTConstMatrix = constMatrix.times(constMatrix.transpose());
    int rank = constMatrix.rank();
    if (rank < initConstraintMatrix.length) { // redundant constraints are specified
      double[][] rMatVal = constMatrixTConstMatrix.qr().getR().getArray();
      List<Double> diag = IntStream.range(0, rMatVal.length).mapToDouble(x->Math.abs(rMatVal[x][x])).boxed().collect(Collectors.toList());
      int[] sortedIndices = IntStream.range(0, diag.size()).boxed().sorted((i, j) -> diag.get(i).compareTo(diag.get(j))).mapToInt(ele->ele).toArray();
      List<Integer> duplicatedEleIndice = IntStream.range(0, diag.size()-rank).map(x -> sortedIndices[x]).boxed().collect(Collectors.toList());
      return genRedundantConstraint(state, duplicatedEleIndice);
    }
    return null;
  }
  
  public static List<String> genRedundantConstraint(ComputationState state, List<Integer> duplicatedEleIndics) {
    List<String> redundantConstraint = new ArrayList<>();
    int numRedEle = duplicatedEleIndics.size();
    for (Integer redIndex : duplicatedEleIndics)
      redundantConstraint.add(grabRedundantConstraintMessage(state, redIndex));

    return redundantConstraint;
  }

  public static String grabRedundantConstraintMessage(ComputationState state, Integer constraintIndex) {
    // figure out which constraint among state._fromBetaConstraints, state._equalityConstraints, 
    // state._lessThanEqualToConstraints is actually redundant
    LinearConstraints redundantConst = getConstraintFromIndex(state, constraintIndex);
    if (redundantConst != null) {
      boolean isBetaConstraint = redundantConst._constraints.size() < 3 ? true : false;
      boolean standardize = state.activeData()._normMul != null ? true : false;
      StringBuilder sb = new StringBuilder();
      DataInfo dinfo = state.activeData();
      List<String> trainNames = Arrays.stream(dinfo._adaptedFrame.names()).collect(Collectors.toList());
      sb.append("This constraint is redundant ");
      
      if (isBetaConstraint && standardize) {  // beta constraint and standardization is needed
        double constantVal = 0;
        int colIndex = 0;
        for (String coefName : redundantConst._constraints.keySet()) {
          double constrVal = redundantConst._constraints.get(coefName);
          if ("constant".equals(coefName)) {
            constantVal = constrVal;
          } else {
            colIndex = trainNames.indexOf(coefName)-dinfo._cats;
            sb.append(coefName);  // coefficient value is always 1.
          }
        }
        sb.append(constantVal*dinfo._normMul[colIndex]);
      } else {  // beta constraint without standardization or linear constraint with or without standardization
        for (String coefName : redundantConst._constraints.keySet()) {
          double constrVal = redundantConst._constraints.get(coefName);
          if (constrVal != 0) {
            if (constrVal > 0)
              sb.append('+');

            // linearconstraint for numerical columns only that needs standardization
            if (standardize && !"constant".equals(coefName) && trainNames.contains(coefName)) {
              int colInd = trainNames.indexOf(coefName) - dinfo._cats;
              sb.append(constrVal/dinfo._normMul[colInd]);
            } else {  // no standardization or standardization but coefName is constant
              sb.append(constrVal);
            }
            if (!"constant".equals(coefName)) {
              sb.append('*');
              sb.append(coefName);
            }
          }
        }
      }
      sb.append(" <= or == 0.");
      sb.append(" Please remove it from your beta/linear constraints.");
      return sb.toString();
    } else {
      return null;
    }
  }
  
  public static LinearConstraints getConstraintFromIndex(ComputationState state, Integer constraintIndex) {
    int constIndexWOffset = constraintIndex;
    if (state._fromBetaConstraints != null) {
      if (constIndexWOffset < state._fromBetaConstraints.length) {
        return state._fromBetaConstraints[constIndexWOffset];
      } else {
        constIndexWOffset -= state._fromBetaConstraints.length;
      }
    }
    
    if (state._equalityConstraints != null) {
      if (constIndexWOffset < state._equalityConstraints.length) {
        return state._equalityConstraints[constIndexWOffset];
      } else {
        constIndexWOffset -= state._equalityConstraints.length;
      }
    }
    
    if (state._lessThanEqualToConstraints != null && constIndexWOffset < state._lessThanEqualToConstraints.length) {
      return state._lessThanEqualToConstraints[constIndexWOffset];
    }
    return null;
  }
}
