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
    extractConstraintCoeffs(state, constraintNamesList);
    // form double matrix
    int constraintNameLen = constraintNamesList.size();
    double[][] initConstraintMatrix = new double[constraintNameLen][constraintNameLen];
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
  
  public static List<String> extractConstraintCoeffs(ComputationState state, List<String> constraintCoeffName) {
    if (state._fromBetaConstraints != null)   // extract coefficients constraints
      extractCoeffNames(constraintCoeffName, state._fromBetaConstraints);

    if (state._equalityConstraints != null)
      extractCoeffNames(constraintCoeffName, state._equalityConstraints);

    if (state._lessThanEqualToConstraints != null)
      extractCoeffNames(constraintCoeffName, state._lessThanEqualToConstraints);
    
    // remove duplicates in the constraints names
    Set<String> noDuplicateNames = new HashSet<>(constraintCoeffName);
    return new ArrayList<>(noDuplicateNames);
  }
  
  public static List<String> extractCoeffNames(List<String> coeffList, LinearConstraints[] constraints) {
    int numConst = constraints.length;
    for (int index=0; index<numConst; index++) {
      Set<String> keys = constraints[index]._constraints.keySet();
      coeffList.addAll(keys);
    }
    return coeffList;
  }
  
  public static boolean foundRedundantConstraints(final double[][] initConstraintMatrix) {
    Matrix constMatrix = new Matrix(initConstraintMatrix);
    return false;
  }
}
