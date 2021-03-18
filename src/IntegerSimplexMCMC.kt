import lib.SettableLazy
import lib.sparseVector.SparseVector
import lib.sparseVector.asVector
import org.apache.commons.math3.util.CombinatoricsUtils
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashSet
import kotlin.math.*
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
// TODO: Pivots tend to reduce the degeneracy probability significantly, causing
//  high rejection ratios and problems coming out of fractional states.
//  Could improve by doing row swaps after each pivot?
open class IntegerSimplexMCMC<T>:
GridMapSimplex<T>
//    Simplex<T,GridMapMatrix<T>>, FieldOperators<T>
//    RowSimplex<T>, FieldOperators<T>
 where T : Comparable<T>, T : Number {

    val logPmf: (SparseVector<T>) -> Double

    // parameters
    val degeneratePivotWeight: Double       = 0.01
    val probOfRowSwap: Double               = 0.01
    val fractionPenaltyK: Double            = -9.0

    // cached values
    val columnPivotLimits =
        ArrayList<T>(M.nCols - 1)      // the upper limit for the value of each column (for easy identification of pivots)
    val columnWeights = MutableCategoricalArray(M.nCols - 1)      // weighted sum of number of pivots in each column
    var state: Cache
    var revertState: PivotReverter<T>
//    var proposedLogFractionalPenalty: Double

    // sample statistics
    var nSamples: Int = 0
    var nFractionalSamples: Int = 0

    inner class Cache {
        var currentSample: SparseVector<T>  by SettableLazy { X(false) }
        var logDegeneracyProb: Double       by SettableLazy { logDegeneracyProb() }
        var logPX: Double                  by SettableLazy { this@IntegerSimplexMCMC.logPmf(currentSample) }
        val logFractionalPenalty: Double
        val logProbOfPivotState: Double
            get() = logPX + logDegeneracyProb + logFractionalPenalty


        constructor(logFractionalPenalty: Double) {
            this.logFractionalPenalty = logFractionalPenalty
        }

        constructor(currentSample: SparseVector<T>, logDegeneracyProb: Double, logPmf: Double, logFractionalPenalty: Double) {
            this.currentSample = currentSample
            this.logDegeneracyProb = logDegeneracyProb
            this.logPX = logPmf
            this.logFractionalPenalty = logFractionalPenalty
        }
    }

    class PivotReverter<T>(val reversePivot: PivotPoint, val cache: IntegerSimplexMCMC<T>.Cache) where T : Comparable<T>, T : Number


    constructor(
        constraints: List<MutableConstraint<T>>,
        initialSample: SparseVector<T>,
        logPmf: (SparseVector<T>) -> Double
    ) :
            super(constraints, emptyMap<Int, T>().asVector(initialSample.operators), initialSample) {
        this.logPmf = logPmf
        for (col in 0 until M.nCols - 1) {
            columnPivotLimits.add(zero)
            updatePivotState(col)
        }
        state = Cache(0.0)
        revertState = PivotReverter(PivotPoint(0,0),state)
    }


    // Choose a pivot
    // and reject based on Metropolis-Hastings
    // returns the next sample
    fun nextSample(): SparseVector<T> {
        if (Random.nextDouble() < probOfRowSwap) {
            processRowSwapProposal(Random.nextInt(M.nRows - 1), Random.nextInt(M.nRows - 1))
        } else {
            var proposalPivot = proposePivot()
            processProposal(proposalPivot)
        }
        ++nSamples
        if(state.logFractionalPenalty != 0.0) {
            ++nFractionalSamples
            println("${nSamples} In fractional state. Fractional penalty ${state.logFractionalPenalty}")
        }
        return state.currentSample
    }


    fun processRowSwapProposal(row1: Int, row2: Int) {
        swapRows(row1,row2)
        val newLogDegeneracyProb = logDegeneracyProb()
        val acceptance = exp(min(0.0, newLogDegeneracyProb - state.logDegeneracyProb))
        if(Random.nextDouble() < acceptance) {
            state = Cache(state.currentSample, newLogDegeneracyProb, state.logPX, state.logFractionalPenalty)
        } else {
            swapRows(row2,row1)
        }
    }


    fun processProposal(proposalPivot: PivotPoint) {
        val acceptanceDenominator = state.logProbOfPivotState + ln(transitionProb(proposalPivot))
        forwardPivot(proposalPivot)
//        forwardPivotWithRowSwap(proposalPivot)
        val acceptanceNumerator = state.logProbOfPivotState + ln(transitionProb(revertState.reversePivot))
        val logAcceptance = min(acceptanceNumerator - acceptanceDenominator, 0.0)
        if (!logAcceptance.isNaN() && Random.nextDouble() >= exp(logAcceptance)) { // explicity accept if both numerator and denominator are -infinity
            if(revertState.cache.logFractionalPenalty != 0.0 && state.logFractionalPenalty == 0.0) {
                println("Rejecting fraction->integer transition with denominator = ${revertState.cache.logPX} + ${revertState.cache.logDegeneracyProb} + ${revertState.cache.logFractionalPenalty} + ${acceptanceDenominator - revertState.cache.logFractionalPenalty - revertState.cache.logDegeneracyProb - revertState.cache.logPX} = $acceptanceDenominator, numerator = ${state.logPX} + ${state.logDegeneracyProb} + ${ln(transitionProb(revertState.reversePivot))} = $acceptanceNumerator, acceptance = $logAcceptance")
            }
            if(revertState.cache.logFractionalPenalty == 0.0 && state.logFractionalPenalty != 0.0) {
                println("Rejecting integer->fraction transition with denominator = ${revertState.cache.logPX} + ${revertState.cache.logDegeneracyProb} + ${acceptanceDenominator - revertState.cache.logFractionalPenalty - revertState.cache.logDegeneracyProb - revertState.cache.logPX} = $acceptanceDenominator, numerator = ${state.logPX} + ${state.logDegeneracyProb} +  + ${revertState.cache.logFractionalPenalty} + ${ln(transitionProb(revertState.reversePivot))} = $acceptanceNumerator, acceptance = $logAcceptance")
            }
            revertLastPivot()
//            revertLastPivotWithRowSwap()
        } else {
            // Accept proposal
        }
    }


    fun nextIntegerSample(): SparseVector<T> {
        var integerSample: SparseVector<T>
        var attempts = 0
        do {
            integerSample = nextSample()
            ++attempts
            if(attempts.rem(10) == 0) println("stuck on fraction $attempts. logFractional penalty ${state.logFractionalPenalty}")
        } while(state.logFractionalPenalty != 0.0)
        if(attempts > 1) println("Made $attempts attempts before getting integer sample")
        return integerSample
    }



    fun <R> expectation(nSamples: Int, initialExpectation: R, expectationAccumulator: (SparseVector<T>, R) -> R): R {
        var e = initialExpectation
        var oldSample: SparseVector<T>? = null
        var rejections = 0
        for (s in 1..nSamples) {
            val newSample = nextSample()
            e = expectationAccumulator(newSample, e)
            if (oldSample === newSample) ++rejections
            oldSample = newSample
        }
        println("Rejection ratio = ${rejections.toDouble() / nSamples}")
        return e
    }


    fun forwardPivot(pivot: PivotPoint, newLogFractionalPenalty: Double = logFractionalPenaltyAfterPivot(pivot.col)) {
        revertState = PivotReverter(PivotPoint(pivot.row, basicColsByRow[pivot.row]), state)
        super.pivot(pivot)
        updatePivotInfo(revertState.reversePivot)
        state = Cache(newLogFractionalPenalty)
    }


    fun revertLastPivot() {
        val thisPivot = PivotPoint(revertState.reversePivot.row, basicColsByRow[revertState.reversePivot.row])
        super.pivot(revertState.reversePivot)
        updatePivotInfo(thisPivot)
        state = revertState.cache
    }


    fun swapRows(row1: Int, row2: Int) {
        if (row1 != row2) { // Do swap
            M.rows[row2].weightedPlusAssign(M.rows[row1], -one)
            M.rows[row1].weightedPlusAssign(M.rows[row2], one)
            M.rows[row2].weightedPlusAssign(M.rows[row1], -one)
            M.rows[row2] *= -one
            val basicCol1 = basicColsByRow[row1]
            basicColsByRow[row1] = basicColsByRow[row2]
            basicColsByRow[row2] = basicCol1
        }
    }


    fun updatePivotInfo(pivotedOut: PivotPoint) {
        val colsToUpdate = HashSet<Int>()
        colsToUpdate.addAll(M.rows[pivotedOut.row].nonZeroEntries.keys)
        if (!B[pivotedOut.row].isZero()) {
            for (i in M.columns[pivotedOut.col].nonZeroEntries.keys) {
                if (i != objectiveRow) colsToUpdate.addAll(M.rows[i].nonZeroEntries.keys)
            }
        }
        colsToUpdate.remove(bColumn)
        for (j in colsToUpdate) {
            updatePivotState(j)
        }
    }


    fun updatePivotState(column: Int) {
        val pivotRows = pivotableRows(column, false)
        columnPivotLimits[column] = if (pivotRows.isEmpty()) zero else B[pivotRows[0]] / M[pivotRows[0], column]
        columnWeights[column] =
            if (columnPivotLimits[column] == zero)
                degeneratePivotWeight * pivotRows.size
            else
                pivotRows.size.toDouble()

    }


    fun proposePivot(): PivotPoint {
        // first choose column with weight given by columnWeights
        val col = columnWeights.sample()
        val nPivot = Random.nextInt(nPivots(col))
        val row = M.columns[col].nonZeroEntries.asSequence()
            .filter { it.key != objectiveRow && (B[it.key] / it.value) as Number == columnPivotLimits[col] }
            .drop(nPivot)
            .first()
            .key
        return PivotPoint(row, col)
    }


    // number of pivot points in given column
    fun nPivots(column: Int): Int {
        return if (columnPivotLimits[column] == zero)
            (columnWeights[column] / degeneratePivotWeight).roundToInt()
        else
            columnWeights[column].roundToInt()
    }


    // The probability of choosing this pivot as a transition
    fun transitionProb(pivot: PivotPoint): Double {
        return (1.0 - probOfRowSwap) * columnWeights.P(pivot.col) / nPivots(pivot.col)
    }


    fun logDegeneracyProb(): Double {
        val possiblePivotCols = HashSet<Int>()
        var logProb = 0.0
        var degeneracy = 0
        for (i in M.nRows - 2 downTo 0) {
            if (B[i].isZero()) {
                possiblePivotCols.addAll(M.rows[i].nonZeroEntries.keys)
                logProb -= ln(possiblePivotCols.size.toDouble())
                ++degeneracy
            }
        }
        return logProb +
                CombinatoricsUtils.factorialLog(basicColsByRow.size - degeneracy) -
                CombinatoricsUtils.factorialLog(basicColsByRow.size)
    }


    fun T.roundingDistance(): Double {
        return this.toDouble().roundToInt() - this.toDouble()
    }


    fun logFractionalPenaltyAfterPivot(column: Int): Double {
        val pivotDelta = columnPivotLimits[column]
        return fractionPenaltyK * (M.columns[column].nonZeroEntries.entries
            .sumByDouble { (B[it.key] - pivotDelta*it.value).roundingDistance().absoluteValue } +
                B.nonZeroEntries.entries
                    .filter { M[it.key, column].isZero() }
                    .sumByDouble { it.value.roundingDistance().absoluteValue })
    }


    /************************ TEST STUFF *********************/


    // Trying row swapping after a pivot to see if this keeps the
    // degeneracy prob more stable.
    // In row order, swap newly non-zero rows with
    // the newly zero rows until no more swaps exist
    fun forwardPivotWithRowSwap(pivot: PivotPoint) {
        val oldBKeys = B.nonZeroEntries.keys.toIntArray()
        forwardPivot(pivot)
        doRowSwap(oldBKeys)
    }

    fun revertLastPivotWithRowSwap() {
        val oldBKeys = B.nonZeroEntries.keys.toIntArray()
        revertLastPivot()
        doRowSwap(oldBKeys)
    }

    fun doRowSwap(oldBKeys: IntArray) {
        val newBKeys = B.nonZeroEntries.keys.toIntArray()
        oldBKeys.sort()
        newBKeys.sort()
        var oldIndex = 0
        var newIndex = 0
        var doSwap: Boolean
        do {
            while(oldIndex < oldBKeys.size && newBKeys.binarySearch(oldBKeys[oldIndex]) >= 0) ++oldIndex
            while(newIndex < newBKeys.size && oldBKeys.binarySearch(newBKeys[newIndex]) >= 0) ++newIndex
            if(oldIndex < oldBKeys.size && newIndex < newBKeys.size) swapRows(oldBKeys[oldIndex], newBKeys[newIndex])
            ++oldIndex
            ++newIndex
        } while(oldIndex < oldBKeys.size && newIndex < newBKeys.size)
    }


    fun fractionOfAffectedColumns(pivot: PivotPoint): Double {
        val affectedCol = HashSet<Int>()
        for (i in M.columns[pivot.col].nonZeroEntries.keys) {
            for (j in M.rows[i].nonZeroEntries.keys) {
                affectedCol.add(j)
            }
        }
        return affectedCol.size * 1.0 / M.nCols
    }

    fun checkConsistency() {
        for (col in 0 until M.nCols - 1) {
            val oldLimit = columnPivotLimits[col]
            val oldWeight = columnWeights[col]
            updatePivotState(col)
//            println("new weight = ${columnWeights[col]} old weight = $oldWeight")
            assert(columnPivotLimits[col] == oldLimit)
            assert((columnWeights[col] - oldWeight).absoluteValue < 1e-9)
        }
    }

    fun bestFractionalPenaltyPivot(): Double {
        return (0 until M.nCols-1)
            .filter { isBasicColumn(it) }
            .map { logFractionalPenaltyAfterPivot(it) }
            .max()!!
    }
}

