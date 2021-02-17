import lib.abstractAlgebra.FieldOperators
import lib.sparseMatrix.SparseMatrix
import lib.vector.SparseVector
import lib.vector.asVector
import org.apache.commons.math3.fraction.Fraction
import kotlin.collections.ArrayList
import kotlin.collections.HashSet
import kotlin.math.absoluteValue
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.roundToInt
import kotlin.random.Random

// Adds MCMC pivoting to a Simplex so that the probability of
// being in a state with solution X is the supplied probability
// mass function pmf(X).
//
// In addition to the solution state, X, we add a degeneracy state
// which consists of a set of M.nRows-|X| columns that are zero in X.
// Each degeneracy state maps to a set of pivoted-in columns for the degenerate
// rows via toPivotState(). This should be an external node on the degeneracy graph.
//
//
//
open class SimplexMCMC<T> : Simplex<T> where T: Comparable<T>, T: Number {

    val logPmf: (SparseVector<T>) -> Double

    val fractionPenaltyK: Double          = -1.0
    val degeneratePivotWeight             = 0.001
    val probOfRowSwap                     = 0.01

    var logProbOfPivotState: Double
    val columnPivotLimits = ArrayList<T>(M.nCols - 1)      // the upper limit for the value of each column (for easy identification of pivots)
    val columnWeights = MutableCategoricalArray(M.nCols - 1)      // weighted sum of number of pivots in each column

    init {
        for (col in 0 until M.nCols - 1) {
            columnPivotLimits.add(zero)
            updatePivotState(col)
        }
    }


    constructor(xCoefficients: SparseMatrix<T>, constants: SparseVector<T>, logPmf: (SparseVector<T>) -> Double) :
            super(xCoefficients, constants, emptyMap<Int, T>().asVector(xCoefficients.operators)) {
        this.logPmf = logPmf
        logProbOfPivotState = calcLogProbOfPivotState(X())
    }

    constructor(operators: FieldOperators<T>, constraints: List<Constraint<T>>, logPmf: (SparseVector<T>) -> Double) :
            super(constraints, emptyMap<Int, T>().asVector(operators)) {
        this.logPmf = logPmf
        logProbOfPivotState = calcLogProbOfPivotState(X())
    }

    fun calcLogProbOfPivotState(x: SparseVector<T>) = logPmf(x) + logDegeneracyProb() + logFractionPenalty(x)


    // Choose a positive pivot with uniform probability
    // and reject based on Metropolis-Hastings
    // returns true if the proposal was accepted
    fun mcmcTransition(): Boolean {
        if(Random.nextDouble() < probOfRowSwap) { // row swap never rejects
            swapRows(Random.nextInt(M.nRows-1), Random.nextInt(M.nRows-1))
            return true
        }
        var proposalPivot = proposePivot()
        println("Proposal is degenerate ${isDegenerate(proposalPivot)}")
        val rejectionPivot = PivotPoint(proposalPivot.row, basicColsByRow[proposalPivot.row])
        val originalLogProbOfPivotState = logProbOfPivotState
        val originalLogProbOfTransition = ln(trasitionProb(proposalPivot))
//        println("Fraction of pivot-affected cols ${fractionOfAffectedColumns(proposalPivot)}")
        pivot(proposalPivot)
        updatePivotStates(rejectionPivot)
        val newX = X()
        logProbOfPivotState  = logPmf(newX) + logDegeneracyProb() + logFractionPenalty(newX)
        val logAcceptance =
            logProbOfPivotState + ln(trasitionProb(rejectionPivot)) - originalLogProbOfPivotState - originalLogProbOfTransition
        println("Acceptance = ${exp(logAcceptance)}")
        return if(Random.nextDouble() >= exp(logAcceptance)) {
            pivot(rejectionPivot)
            updatePivotStates(proposalPivot)
            logProbOfPivotState = originalLogProbOfPivotState
            false
        } else {
            true
        }
    }


    fun swapRows(row1: Int, row2: Int) {
        if(row1 != row2) { // Do swap
            M.rows[row2].weightedPlusAssign(M.rows[row1], -one)
            M.rows[row1].weightedPlusAssign(M.rows[row2], one)
            M.rows[row2].weightedPlusAssign(M.rows[row1], -one)
            M.rows[row2] *= -one
        }
    }


    fun updatePivotStates(pivotedOut: PivotPoint) {
        val colsToUpdate = HashSet<Int>()
        colsToUpdate.addAll(M.rows[pivotedOut.row].nonZeroEntries.keys)
        if(!B[pivotedOut.row].isZero()) {
            for(i in M.columns[pivotedOut.col].nonZeroEntries.keys) {
                if(i != objectiveRow) colsToUpdate.addAll(M.rows[i].nonZeroEntries.keys)
            }
        }
        colsToUpdate.remove(bColumn)
        for(j in colsToUpdate) {
            updatePivotState(j)
        }
    }


    fun updatePivotState(column: Int) {
        val pivotRows = pivotableRows(column)
        columnPivotLimits[column] = if(pivotRows.isEmpty()) zero else B[pivotRows[0]]/M[pivotRows[0],column]
        columnWeights[column] =
            if(columnPivotLimits[column] == zero)
                degeneratePivotWeight*pivotRows.size
            else
                pivotRows.size.toDouble()

    }


    fun proposePivot(): PivotPoint {
        // first choose column with weight given by columnWeights
        val col = columnWeights.sample()
        val nPivot = Random.nextInt(nPivots(col))
        val row = M.columns[col].nonZeroEntries.asSequence()
            .filter { it.key != objectiveRow && (B[it.key]/it.value) as Number == columnPivotLimits[col] }
            .drop(nPivot)
            .first()
            .key
        return PivotPoint(row,col)
    }


    // number of pivot points in given column
    fun nPivots(column: Int): Int {
        return if(columnPivotLimits[column] == zero)
            (columnWeights[column]/degeneratePivotWeight).roundToInt()
        else
            columnWeights[column].roundToInt()
    }


    // The probability of choosing this pivot as a transition
    fun trasitionProb(pivot: PivotPoint): Double {
        return (1.0 - probOfRowSwap) * columnWeights.P(pivot.col) / nPivots(pivot.col)
    }


    fun logDegeneracyProb(): Double {
        val possiblePivotCols = HashSet<Int>()
        var logProb = 0.0
        for(i in M.nRows-2 downTo 0) {
            if(B[i].isZero()) {
                possiblePivotCols.addAll(M.rows[i].nonZeroEntries.keys)
                logProb += ln(possiblePivotCols.size.toDouble())
            }
        }
        return logProb
    }


    // The probability of a pivot state is multiplied by this amount if the
    // solution is not on the integer grid.
    //
    // Returns the L1 distance between this point and its rounding, multiplied
    // by fractionPenaltyK
    fun logFractionPenalty(x: SparseVector<T>): Double {
        return fractionPenaltyK * x.nonZeroEntries.values.sumByDouble {
            val xi = it.toDouble()
            (xi - xi.roundToInt()).absoluteValue
        }
    }


    /************************ TEST STUFF *********************/

    fun fractionOfAffectedColumns(pivot: PivotPoint): Double {
        val affectedCol = HashSet<Int>()
        for(i in M.columns[pivot.col].nonZeroEntries.keys) {
            for(j in M.rows[i].nonZeroEntries.keys) {
                affectedCol.add(j)
            }
        }
        return affectedCol.size*1.0/M.nCols
    }

    fun checkConsistency() {
        for(col in 0 until M.nCols-1) {
            val oldLimit = columnPivotLimits[col]
            val oldWeight = columnWeights[col]
            updatePivotState(col)
//            println("new weight = ${columnWeights[col]} old weight = $oldWeight")
            assert(columnPivotLimits[col] == oldLimit)
            assert((columnWeights[col] - oldWeight).absoluteValue < 1e-9)
        }
    }

}