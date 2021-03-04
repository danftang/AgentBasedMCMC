import lib.abstractAlgebra.FractionOperators
import lib.isInteger
import lib.unaryMinus
import lib.plus
import lib.minus
import lib.sparseVector.SparseVector
import lib.sparseVector.asVector
import org.apache.commons.math3.fraction.Fraction
import java.lang.Integer.max
import kotlin.math.absoluteValue
import kotlin.math.roundToInt

class ABMCMC<AGENT : Agent<AGENT>, ACT : Ordered<ACT>> {

    val fractionPenaltyK: Double          = -4.0

    val simplex: SimplexMCMC<Fraction>
    val model: ABM<AGENT, ACT>
    val observations: List<Observation<AGENT, ACT>>


    constructor(
        model: ABM<AGENT, ACT>,
        nTimesteps: Int,
        observations: List<Observation<AGENT, ACT>>,
        initialSample: SparseVector<Fraction>? = null
    ) {
        this.model = model
        this.observations = observations
        val constraints = constraints(model, nTimesteps, observations)
        this.simplex = SimplexMCMC(
            constraints,
            initialSample
                ?: run {
                    ORTools.IntegerSolve(constraints)
                        .withIndex()
                        .filter { it.value != 0.0 }
                        .associate { Pair(it.index, Fraction(it.value)) }
                        .asVector(FractionOperators)
                },
            this::logProb
        )
        assert(simplex.isFullyPivoted())
        assert(simplex.isPrimalFeasible())
    }


    fun logProb(X: SparseVector<Fraction>): Double {
        val trajectory = Trajectory(model, X)
        val prior = trajectory.logPrior()
        val likelihood = observations.sumByDouble { it.logLikelihood(trajectory) }
        val logP = prior + likelihood + logFractionPenalty(X)
//        println("got logprob $prior + $likelihood = $logP")
        return logP
    }

    // The probability of a pivot state is multiplied by this amount if the
    // solution is not on the integer grid.
    //
    // Returns the L1 distance between this point and its rounding, multiplied
    // by fractionPenaltyK
    fun logFractionPenalty(x: SparseVector<Fraction>): Double {
        return fractionPenaltyK * x.nonZeroEntries.values.sumByDouble {
            val xi = it.toDouble()
            (xi - xi.roundToInt()).absoluteValue
        }
    }


    fun nextSample(): Trajectory<AGENT,ACT> {
        return Trajectory(model, nextSampleVector())
    }


    fun nextSampleVector(): SparseVector<Fraction> {
        var integerSample: SparseVector<Fraction>
        var attempts = 0
        do {
            integerSample = simplex.nextSample()
            ++attempts
            if(attempts.rem(10) == 0) println("stuck on fraction $attempts. Penalty is ${logFractionPenalty(integerSample)}")
        } while(!integerSample.isInteger())
        if(attempts > 1) println("Made $attempts attempts before getting integer sample")
        return integerSample
    }


    fun<R> expectation(nSamples: Int, initialExpectation: R, expectationAccumulator: (SparseVector<Fraction>, R) -> R): R {
        var e = initialExpectation
        var oldSample: SparseVector<Fraction>? = null
        var rejections = 0
        for(s in 1..nSamples) {
            val newSample = nextSampleVector()
            e = expectationAccumulator(newSample, e)
            if(oldSample === newSample) ++rejections
            oldSample = newSample
            if(s.rem(10) == 0) println("Got sample $s, ${largestDenominator()}, ${simplex.M.sparsity()}")
        }
        println("Rejection ratio = ${rejections.toDouble()/nSamples}")
        return e
    }


    fun largestDenominator(): Int {
        var maxDenom = 1
        for(entry in simplex.M.nonZeroEntries) {
            maxDenom = max(entry.value.denominator, maxDenom)
            assert(entry.value != Fraction.ZERO)
        }
        return maxDenom
    }

    companion object {
        fun <AGENT : Agent<AGENT>, ACT : Ordered<ACT>> constraints(
            model: ABM<AGENT, ACT>,
            nTimesteps: Int,
            observations: List<Observation<AGENT, ACT>>
        ): List<Constraint<Fraction>> {
            return validTrajectoryConstraints(model, nTimesteps, true) + observations.flatMap { it.eventConstraints() }
        }

        fun <AGENT : Agent<AGENT>, ACT : Ordered<ACT>> validTrajectoryConstraints(
            model: ABM<AGENT, ACT>,
            nTimesteps: Int,
            fermionic: Boolean
        ): List<Constraint<Fraction>> {
            println("Constructing model constraints...")
            val constraints = ArrayList<Constraint<Fraction>>(model.actDomain.size * model.agentDomain.size * 2)
            constraints.addAll(continuityConstraints(nTimesteps, model))
            if(fermionic) constraints.addAll(fermionicConstraints(nTimesteps, model))
            for (state in 0 until model.agentDomain.size) {
                for (act in 0 until model.actDomain.size) {
                    for (t in 0 until nTimesteps) {
                        val event = ABM.Event(t, model.agentDomain[state], model.actDomain[act])
                        constraints.addAll(fermionicXImpliesY(event.ordinal, model.timestepEventConstraints(event)))
                    }
                }
            }
//            println("Valid trajectory constraints are $constraints")
            println("Done")
            return constraints
        }

        fun <AGENT : Agent<AGENT>, ACT : Ordered<ACT>> continuityConstraints(
            nTimesteps: Int,
            abm: ABM<AGENT, ACT>
        ): List<Constraint<Fraction>> {
            val nStates = abm.agentDomain.size
//        val acts = abm.actDomain
            val nActs = abm.actDomain.size
            val constraints = ArrayList<Constraint<Fraction>>((nTimesteps - 1) * nStates)

            // first do leaving edges
            for (t in 1 until nTimesteps) {
                for (state in 0 until nStates) {
                    val coeffs = HashMap<Int, Fraction>()
                    val coeffBase = (t * nStates + state) * nActs
                    for (act in 0 until nActs) {
                        coeffs[coeffBase + act] = Fraction(-1)
                    }
                    constraints.add(Constraint(coeffs, "==", Fraction.ZERO))
                }
            }
            // now do incoming edges
            for (state in 0 until nStates) {
                for (act in 0 until nActs) {
                    val consequences = abm.consequences(abm.agentDomain[state], abm.actDomain[act])
                    for ((resultState, n) in consequences.entries) {
                        val resultIndex = resultState.ordinal
                        for (t in 0 until nTimesteps - 1) {
                            constraints[t * nStates + resultIndex]
                                .coefficients[(t * nStates + state) * nActs + act] = Fraction(n)
                        }
                    }
                }
            }
//            println("Continuity constraints are $constraints")
            return constraints
        }


        fun <AGENT : Agent<AGENT>, ACT : Ordered<ACT>> fermionicConstraints(
            nTimesteps: Int,
            abm: ABM<AGENT, ACT>
        ): List<Constraint<Fraction>> {
            val nStates = abm.agentDomain.size
            val nActs = abm.actDomain.size
            val constraints = ArrayList<Constraint<Fraction>>((nTimesteps + 1) * nStates)

            for (t in 0 until nTimesteps) {
                for (state in 0 until nStates) {
                    val coeffs = HashMap<Int, Fraction>()
                    val coeffBase = (t * nStates + state) * nActs
                    for (act in 0 until nActs) {
                        coeffs[coeffBase + act] = Fraction(1)
                    }
                    constraints.add(Constraint(coeffs, "<=", Fraction.ONE))
                }
            }

            val coeffsByState = Array(nStates) { HashMap<Int, Fraction>() }
            for (state in 0 until nStates) {
                val coeffBase = ((nTimesteps - 1) * nStates + state) * nActs
                for (act in 0 until nActs) {
                    val consequences = abm.consequences(abm.agentDomain[state], abm.actDomain[act])
                    for ((resultState, n) in consequences.entries) {
                        coeffsByState[resultState.ordinal][coeffBase + act] = Fraction(n)
                    }
                }
            }
            for (state in 0 until nStates) {
                constraints.add(Constraint(coeffsByState[state], "<=", Fraction.ONE))
            }
//            println("Fermionic constraints are $constraints")
            return constraints
        }

        // Returns constraint x -> y
        // under the assumption that
        // 0 <= x <= 1
        // and 0 <= y_i <= 1
        // by using the identity
        //
        fun fermionicXImpliesY(x: Int, y: Constraint<Fraction>): List<Constraint<Fraction>> {
            return if (y.relation == "==") {
                fermionicXImpliesY(x, Constraint(y.coefficients, "<=", y.constant)) +
                        fermionicXImpliesY(x, Constraint(y.coefficients, ">=", y.constant))

            } else {
                val const: Fraction
                val coeffs = if (y.relation == "<=") {
                    const = y.constant
                    HashMap(y.coefficients)
                } else {
                    const = -y.constant
                    val negatedCoeffs = HashMap<Int, Fraction>(y.coefficients.size)
                    y.coefficients.mapValuesTo(negatedCoeffs) { -it.value }
                    negatedCoeffs
                }
                var maxVal = Fraction.ZERO
                coeffs.values.forEach { if (it > Fraction.ZERO) maxVal += it }
                coeffs[x] = maxVal - const
                listOf(Constraint(coeffs, "<=", maxVal))
            }
        }

        fun fermionicXImpliesY(x: Int, y: List<Constraint<Fraction>>): List<Constraint<Fraction>> =
            y.flatMap { fermionicXImpliesY(x, it) }


        // converts a constraint in terms of state occupation numbers into a constraint on acts
        // in a given timestep
        fun <AGENT : Agent<AGENT>, ACT : Ordered<ACT>> stateConstraintToActConstraint(
            stateConstraint: Constraint<Fraction>,
            timestep: Int,
            abm: ABM<AGENT, ACT>
        ): Constraint<Fraction> {
            val actCoeffs = HashMap<Int, Fraction>()
            val nActs = abm.actDomain.size
            val timestepBase = abm.agentDomain.size * nActs * timestep
            stateConstraint.coefficients.forEach { (state, coefficient) ->
                for (act in 0 until nActs) {
                    actCoeffs[timestepBase + state * nActs + act] = coefficient
                }
            }
            return Constraint(actCoeffs, stateConstraint.relation, stateConstraint.constant)
        }

        fun SparseVector<Fraction>.isInteger(): Boolean {
            return this.nonZeroEntries.all { it.value.isInteger() }
        }


    }
}