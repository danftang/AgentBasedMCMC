import lib.sparseMatrix.SparseMatrix
import lib.vector.SparseVector
import lib.vector.asMapVector
import org.apache.commons.math3.distribution.BinomialDistribution
import java.util.*
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
class SimplexMCMC<T>(
    xCoefficients: SparseMatrix<T>,
    constants: SparseVector<T>,
    val logPmf: (SparseVector<T>) -> Double
) : Simplex<T>(xCoefficients, constants, emptyMap<Int,T>().asMapVector(xCoefficients.operators)) where T: Comparable<T>, T: Number {

    val fractionPenaltyK: Double          = ln(0.5)
    val degeneratePivotWeight             = 0.001

    val binomial: BinomialDistribution    = BinomialDistribution(basicColsByRow.size, 0.5)
//    var positivePivots: List<PivotPoint>  = emptyList()
    var logProbOfX: Double
//    var nDegeneratePivots: Int            = 0
    val columnPivotLimits: ArrayList<T>
    val columnWeights: MutableCategoricalArray
    var degeneracyState: MutableList<Int>

    class PivotState<T>(val limit: T, val nPivots: Int)

    init {
//        updatePivots()
//        positivePivots = allPositivePivotPoints()
        logProbOfX = logProbOf(X())
        columnPivotLimits = ArrayList(M.nCols-1)
        columnWeights = MutableCategoricalArray(M.nCols-1)
        for(col in 0 until M.nCols-1) {
            columnPivotLimits.add(zero)
            updatePivotState(col)
        }
//        degeneracyState = MutableCategoricalMap(MutableCategoricalMap.MapType.TREEMAP)
        degeneracyState = basicColsByRow
            .withIndex()
            .filter { (row, col) -> B[row].isZero() }
            .map { it.value }
            .toMutableList()
    }

    // pivots: For each column, the pivot limit for that column
    // this can be updated incrementally with the pivot
    // choose a column based on the pivot limit
//    inner class PivotLimits {
//        val columnLimits: ArrayList<T?> = ArrayList(M.nCols-1)
//
//        init {
//            for(j in 0 until M.nCols-1) columnLimits.add(columnPivotLimit(j))
//        }
//    }


    // Choose a positive pivot with uniform probability
    // and reject based on Metropolis-Hastings
    // returns true if the proposal was accepted
    fun mcmcTransition(): Boolean {

        var proposalPivot = proposePivot()
        println("Proposal is degenerate ${isDegenerate(proposalPivot)}")
        val rejectionPivot = PivotPoint(proposalPivot.row, basicColsByRow[proposalPivot.row])
        val originalLogProbOfX = logProbOfX
        println("Fraction of pivot-affected cols ${fractionOfAffectedColumns(proposalPivot)}")
        pivot(proposalPivot)
        updatePivotStates(rejectionPivot)
        logProbOfX  = logProbOf(X())
        val logAcceptance = logProbOfX - originalLogProbOfX
        println("Acceptance = ${exp(logAcceptance)}")
        return if(Random.nextDouble() >= exp(logAcceptance)) {
            pivot(rejectionPivot)
            updatePivotStates(proposalPivot)
            logProbOfX = originalLogProbOfX
            false
        } else {
            true
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

    fun nPivots(column: Int): Int {
        return if(columnPivotLimits[column] == zero)
            (columnWeights[column]/degeneratePivotWeight).roundToInt()
        else
            columnWeights[column].roundToInt()
    }


    fun perturbDegenerateState() {
        val stateToDrop = degeneracyState.random()
        val stateToAdd =
    }

//    fun pivotAndUpdateLimits(pivot: PivotPoint) {
//        val modifiedRows = M.columns[pivot.col].nonZeroEntries.asSequence()
//            .map{ it.key }
//            .filter { it != pivot.row }
//            .toList()
//
//        pivot(pivot)
//
//        modifiedRows.forEach { i ->
//            M.rows[i].nonZeroEntries.forEach { (j, Mij) ->
//                if(Mij > zero && B[i]/Mij)
//            }
//        }
//
//        val i = pivot.row
//        val j = pivot.col
//        assert(i>=0 && i < M.nRows-1)
//        assert(j>=0 && j < M.nCols-1)
//        var Mij = M[i,j]
//
//        M.rows[i] /= Mij
//        val colEntries = M.columns[j].nonZeroEntries.asSequence().filter { it.key != i }.toList()
//        val rowEntries = M.rows[i].nonZeroEntries.entries.toList()
//        for(rowEntry in rowEntries) {
//            for(colEntry in colEntries) {
//                val outerProd = rowEntry.value * colEntry.value
//                M.mapAssign(colEntry.key, rowEntry.key) { oldVal -> oldVal - outerProd }
//            }
//        }
//        basicColsByRow[i] = j
//
//    }

//    fun logProbOfProposal(pivot: PivotPoint): Double {
//        val nNonDegeneratePivots = positivePivots.size - nDegeneratePivots
////        println("nonDegeneratePivots = ${nNonDegeneratePivots}/${positivePivots.size}")
//        val pivotWeight = if(isDegenerate(pivot)) degeneratePivotWeight else 1.0
//        return ln(pivotWeight/(nNonDegeneratePivots + nDegeneratePivots*degeneratePivotWeight))
//    }

    fun logFractionPenalty(x: SparseVector<T>): Double {
        return fractionPenaltyK * x.nonZeroEntries.values.count { it.toDouble() != it.toInt().toDouble() }

    }

    fun logProbOf(x: SparseVector<T>): Double =
        logPmf(x) + logFractionPenalty(x) - binomial.logProbability(x.nonZeroEntries.size)

//    fun updatePivots() {
//        positivePivots = allPositivePivotPoints()
//        nDegeneratePivots = positivePivots.count { isDegenerate(it) }
//    }

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