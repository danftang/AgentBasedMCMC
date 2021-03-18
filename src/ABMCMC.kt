import lib.abstractAlgebra.FractionOperators
import lib.isInteger
import lib.unaryMinus
import lib.plus
import lib.minus
import lib.sparseVector.SparseVector
import lib.sparseVector.asVector
import org.apache.commons.math3.fraction.Fraction
import java.lang.Integer.max
import java.time.Instant
import kotlin.math.ln

class ABMCMC<AGENT : Agent<AGENT>, ACT : Ordered<ACT>> {


    val simplex: IntegerSimplexMCMC<Fraction>
    val model: ABM<AGENT, ACT>
    val observations: List<Observation<AGENT, ACT>>


    //////////////// Fractional vertex tests

    var nSamples = 0
    var nRejections = 0
    var fractionalRunLength = 0



    //////////////// End of fractional vertex tests

    constructor(
        model: ABM<AGENT, ACT>,
        nTimesteps: Int,
        observations: List<Observation<AGENT, ACT>>,
        initialSample: SparseVector<Fraction>? = null
    ) {
        this.model = model
        this.observations = observations
        val constraints = constraints(model, nTimesteps, observations)
        this.simplex = IntegerSimplexMCMC(
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

//        currentFractionalPenalty = logFractionPenalty(simplex.X(true))
    }


    fun logProb(X: SparseVector<Fraction>): Double {
        val trajectory = Trajectory(model, X)
        val prior = trajectory.logPrior()
        val likelihood = observations.sumByDouble { it.logLikelihood(trajectory) }
//        val penalty = logFractionPenalty(X)
        val logP = prior + likelihood //+ penalty
//        if(penalty < 0.0) println("Got fractional penalty $penalty")
//        println("prior logprob $prior likelihood logprob $likelihood fraction penalty $penalty = $logP")
        return logP
    }



    fun nextSample(): Trajectory<AGENT,ACT> {
        return Trajectory(model, simplex.nextIntegerSample())
    }




    fun<R> expectation(nSamples: Int, initialExpectation: R, expectationAccumulator: (SparseVector<Fraction>, R) -> R): R {
        var e = initialExpectation
        var oldSample: SparseVector<Fraction>? = null
        var rejections = 0
        var lastTime = Instant.now().toEpochMilli()
        for(s in 1..nSamples) {
            val newSample = simplex.nextIntegerSample()
            e = expectationAccumulator(newSample, e)
            if(oldSample === newSample) ++rejections
            oldSample = newSample
            if(s.rem(100) == 0) {
                val now = Instant.now().toEpochMilli()
                println("Got to sample $s in ${(now-lastTime)/1000.0}s largest Numerator,Denominator ${largestNumeratorDenominator()}, Size ${simplex.entryMap.size}, Degeneracy ${simplex.degeneracy()} logPiv = ${simplex.state.logProbOfPivotState} logPX = ${simplex.state.logPX}, logPDegeneracy = ${simplex.state.logDegeneracyProb}")
                lastTime = now
//                println(simplex.fractionalLogP - simplex.state.logPX - ln(simplex.transitionProb(simplex.proposePivot())))
            }
        }
        println("Rejection ratio = ${rejections.toDouble()/nSamples}")
        return e
    }


    fun largestDenominator(): Int {
        var maxDenom = 1
        for(entry in simplex.entryMap.entries) {
            maxDenom = max(entry.value.denominator, maxDenom)
            assert(entry.value != Fraction.ZERO)
        }
        return maxDenom
    }

    fun largestNumeratorDenominator(): Pair<Int,Int> {
        var maxDenom = 1
        var maxNum = 1
        for(entry in simplex.entryMap.entries) {
            maxDenom = max(entry.value.denominator, maxDenom)
            maxNum = max(entry.value.numerator, maxDenom)
            assert(entry.value != Fraction.ZERO)
        }
        return Pair(maxNum, maxDenom)
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