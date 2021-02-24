import org.apache.commons.math3.fraction.Fraction

interface Observation<AGENT : Agent<AGENT>,ACT : Ordered<ACT>> {
    fun logLikelihood(trajectory: Trajectory<AGENT,ACT>): Double
    fun constraints(): List<Constraint<Fraction>>
}