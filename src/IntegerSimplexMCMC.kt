import lib.SettableLazy
import lib.sparseVector.SparseVector
import lib.sparseVector.asVector
import org.apache.commons.math3.util.CombinatoricsUtils
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
// TODO: Need to ensure fractional->fractional transitions don't wander
// TODO: to very unlikely states. Also, improve mixing by doing row swaps after all pivots?
open class IntegerSimplexMCMC<T> : Simplex<T> where T : Comparable<T>, T : Number {

    val logPmf: (SparseVector<T>) -> Double

    // parameters
    val degeneratePivotWeight: Double       = 0.01
    val probOfRowSwap: Double               = 0.01
    val targetFractionalProportion: Double  = 0.01
    val fractionalPUpdateAlpha: Double      = 0.1
    val fractionPenaltyK: Double            = -1.0
    var fractionalLogP: Double

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
            get() = logPX + logDegeneracyProb


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
        constraints: List<Constraint<T>>,
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
        fractionalLogP = state.logProbOfPivotState + ln(trasitionProb(proposePivot())) // first guess at fractional weight
        println("initialising fractionalLogP = $fractionalLogP")
    }


//    fun calcLogProbOfPivotState(x: SparseVector<T>): Double {
//        val pmf = logPmf(x)
//        val degeneracyP = logDegeneracyProb()
////        println("logPmf $pmf degeneracy $degeneracyP = ${pmf+degeneracyP}")
//        return pmf + degeneracyP
//    }

    // Choose a pivot
    // and reject based on Metropolis-Hastings
    // returns the next sample
    fun nextSample(): SparseVector<T> {
        if (Random.nextDouble() < probOfRowSwap) {
            processRowSwapProposal(Random.nextInt(M.nRows - 1), Random.nextInt(M.nRows - 1))
        } else {
            var proposalPivot = proposePivot()
            val proposedLogFractionalPenalty = logFractionalPenaltyAfterPivot(proposalPivot.col)
            if (proposedLogFractionalPenalty == 0.0) {
                processIntegerProposal(proposalPivot, proposedLogFractionalPenalty)
            } else {
                processFractionalProposal(proposalPivot, proposedLogFractionalPenalty)
            }
        }
        ++nSamples
        if(state.logFractionalPenalty > 0.0) ++nFractionalSamples
//        if(nSamples.rem(1024) == 0) updateFractionalP()
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


    // destination is integer
    fun processIntegerProposal(proposalPivot: PivotPoint, proposedLogFractionalPenalty: Double) {
//        if (isDegenerate(proposalPivot)) println("Proposal is degenerate")
        val acceptanceDenominator = calcLogAcceptanceDenominator(proposalPivot)
        forwardPivot(proposalPivot, proposedLogFractionalPenalty)
        val acceptanceNumerator = state.logProbOfPivotState + ln(trasitionProb(revertState.reversePivot))
        val logAcceptance = min(acceptanceNumerator - acceptanceDenominator, 0.0)
//        println("Log acceptance = $logProbOfPivotState + ${ln(trasitionProb(rejectionPivot))} - $originalLogProbOfPivotState - $originalLogProbOfTransition")
//        println("Acceptance = ${exp(logAcceptance)}")
        // explicity accept if both numerator and denominator are -infinity
        if (!logAcceptance.isNaN() && Random.nextDouble() >= exp(logAcceptance)) {
//            println("Rejecting. pivot state prob ratio ${exp(logPivotStateRatio)} transition ratio ${exp(logTransitionProbRatio)}")
//            println()
            if(revertState.cache.logFractionalPenalty != 0.0 && state.logFractionalPenalty == 0.0) println("Rejecting fraction->integer transition with logAcceptance ${state.logPX}, ${state.logDegeneracyProb}, $acceptanceNumerator - $acceptanceDenominator = $logAcceptance")
            revertLastPivot()
        } else {
//            println("Accepting. pivot state log prob ${logProbOfPivotState}")
        }
//        fractionalLogP = (255.0 * fractionalLogP + state.logProbOfPivotState)/256.0
    }


    // destination is fractional
    fun processFractionalProposal(proposalPivot: PivotPoint, proposedLogFractionalPenalty: Double) {
//        println("In fractional state. Best pivot is ${bestFractionalPenaltyPivot()}")
        val logAcceptance = min(0.0,fractionalLogP + proposedLogFractionalPenalty - calcLogAcceptanceDenominator(proposalPivot))
        if(Random.nextDouble() < exp(logAcceptance)) {
//            println("Accepted fractional proposal with fractional penalty $proposedLogFractionalPenalty")
            forwardPivot(proposalPivot, proposedLogFractionalPenalty)
        } else {
            println("Rejected proposal with fractional penalty $proposedLogFractionalPenalty")
            if(proposedLogFractionalPenalty == 0.0) println("Rejected integer proposal from fractional state")
        }
    }


    fun updateFractionalP() {
        val observedFractionalProportion = if(nFractionalSamples != 0) nFractionalSamples.toDouble()/nSamples else 1.0/nSamples
        fractionalLogP += ln(
            (1.0-fractionalPUpdateAlpha) +
                    fractionalPUpdateAlpha * targetFractionalProportion/observedFractionalProportion
        )
        println("Updating fractional weight. Proportion of fractional samples = ${nFractionalSamples.toDouble()/nSamples}")
        println("New log weight = ${fractionalLogP}")
        nSamples = 0
        nFractionalSamples = 0
    }


    fun calcLogAcceptanceDenominator(proposalPivot: PivotPoint): Double =
        if(state.logFractionalPenalty == 0.0) {
            state.logProbOfPivotState + ln(trasitionProb(proposalPivot))
        } else {
            fractionalLogP + state.logFractionalPenalty
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
    fun trasitionProb(pivot: PivotPoint): Double {
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

    // The probability of a pivot state is multiplied by this amount if the
    // solution is not on the integer grid.
    //
    // Returns the L1 distance between this point and its rounding, multiplied
    // by fractionPenaltyK
//    fun logFractionProb(x: SparseVector<T>): Double {
//        return fractionPenaltyK * x.nonZeroEntries.values.sumByDouble {
//            val xi = it.toDouble()
//            (xi - xi.roundToInt()).absoluteValue
//        } + logMeanPBeta
//    }


//    fun SparseVector<T>.isInteger(): Boolean {
//        return this.nonZeroEntries.all { it.value.isInteger() }
//    }

//    fun T.isInteger(): Boolean {
//        return toDouble() == toInt().toDouble()
//    }

    fun T.roundingDistance(): Double {
        return this.roundToInt() - this.toDouble()
    }


    // returns true if pivot will result in an integer point
    // i.e. if B - M[.,j]B[i]/M[i,j] is integer
//    fun isIntegerPivot(pivot: PivotPoint): Boolean {
//        val Mij = M[pivot.row,pivot.col]
//        if(Mij as Number == one && isInIntegerState) return true // most cases
//        val pivotDelta = B[pivot.row]/Mij
//        return M.columns[pivot.col].nonZeroEntries.all { (B[it.key] - pivotDelta*it.value).isInteger() } &&
//                B.nonZeroEntries.all { it.value.isInteger() || !M[it.key,pivot.col].isZero() }
//    }

    fun logFractionalPenaltyAfterPivot(column: Int): Double {
        val pivotDelta = columnPivotLimits[column]
        return fractionPenaltyK * (M.columns[column].nonZeroEntries.entries
            .sumByDouble { (B[it.key] - pivotDelta*it.value).roundingDistance().absoluteValue } +
                B.nonZeroEntries.entries
                    .filter { M[it.key, column].isZero() }
                    .sumByDouble { it.value.roundingDistance().absoluteValue })
    }


    /************************ TEST STUFF *********************/

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
        return variableIndices
            .filter { isBasicColumn(it) }
            .map { logFractionalPenaltyAfterPivot(it) }
            .max()!!
    }
}

