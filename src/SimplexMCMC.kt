import lib.sparseMatrix.SparseMatrix
import lib.vector.SparseVector
import lib.vector.asMapVector
import org.apache.commons.math3.distribution.BinomialDistribution
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.min
import kotlin.math.pow
import kotlin.random.Random

// Adds MCMC pivoting to a Simplex so that the probability of
// being in a state with solution X is the supplied probability
// mass function pmf(X).
//
//
class SimplexMCMC<T>(
    xCoefficients: SparseMatrix<T>,
    constants: SparseVector<T>,
    val logPmf: (SparseVector<T>) -> Double
) : Simplex<T>(xCoefficients, constants, emptyMap<Int,T>().asMapVector(xCoefficients.operators)) where T: Comparable<T>, T: Number {

    val fractionPenaltyK: Double          = ln(0.5)
    val binomial: BinomialDistribution    = BinomialDistribution(basicColsByRow.size, 0.5)
    var positivePivots: List<PivotPoint>
    var logProbOfX: Double

    init {
        positivePivots = allPositivePivotPoints()
        logProbOfX = logProbOf(X())
    }

    // Choose a positive pivot with uniform probability
    // and reject based on Metropolis-Hastings
    fun mcmcTransition() {
        val proposalPivot = positivePivots.random()
        var logAcceptance = ln(positivePivots.size.toDouble()) - logProbOfX
        val rejectionPivot = PivotPoint(proposalPivot.row, basicColsByRow[proposalPivot.row])
        pivot(proposalPivot)
        val newPositivePivots = allPositivePivotPoints()
        val newLogProbOfX  = logProbOf(X())
        logAcceptance += newLogProbOfX - ln(newPositivePivots.size.toDouble())
        val acceptanceProb = min(1.0, exp(logAcceptance))
        if(Random.nextDouble() >= acceptanceProb) {
            pivot(rejectionPivot)
        } else {
            positivePivots = newPositivePivots
            logProbOfX = newLogProbOfX
        }
    }

    fun logFractionPenalty(x: SparseVector<T>): Double {
        return fractionPenaltyK * x.nonZeroEntries.values.count { it.toDouble() != it.toInt().toDouble() }

    }

    fun logProbOf(x: SparseVector<T>): Double =
        logPmf(x) + logFractionPenalty(x) - binomial.logProbability(x.nonZeroEntries.size)

}