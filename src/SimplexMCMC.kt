import lib.sparseMatrix.SparseMatrix
import lib.vector.SparseVector
import lib.vector.asMapVector
import org.apache.commons.math3.distribution.BinomialDistribution
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
    val pmf: (SparseVector<T>) -> Double
) : Simplex<T>(xCoefficients, constants, emptyMap<Int,T>().asMapVector(xCoefficients.operators)) where T: Comparable<T>, T: Number {

    val fractionPenaltyK: Double          = 0.5
    val binomial: BinomialDistribution    = BinomialDistribution(basicColsByRow.size, 0.5)
    var positivePivots: List<PivotPoint>
    var probabilityOfX: Double

    init {
        positivePivots = allPositivePivotPoints()
        probabilityOfX = probabilityOf(X())
    }

    // Choose a pivot with a probability
    fun mcmcTransition() {
        val proposalPivot = positivePivots.random()
        val acceptanceDenominator = probabilityOfX / positivePivots.size
        val rejectionPivot = PivotPoint(proposalPivot.row, basicColsByRow[proposalPivot.row])
        pivot(proposalPivot)
        val newPositivePivots = allPositivePivotPoints()
        val newProbOfX  = probabilityOf(X())
        val acceptanceNumerator = newProbOfX / newPositivePivots.size
        val acceptanceProb = min(1.0, acceptanceNumerator/acceptanceDenominator)
        if(Random.nextDouble() >= acceptanceProb) {
            pivot(rejectionPivot)
        } else {
            positivePivots = newPositivePivots
            probabilityOfX = newProbOfX
        }
    }

    fun fractionPenalty(x: SparseVector<T>): Double {
        return fractionPenaltyK.pow(
            x.nonZeroEntries.values.count { it.toDouble() != it.toInt().toDouble() }
        )
    }

    fun probabilityOf(x: SparseVector<T>): Double =
        pmf(x) * fractionPenalty(x) / binomial.probability(x.nonZeroEntries.size)

}