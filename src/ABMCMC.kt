import lib.abstractAlgebra.FractionOperators
import lib.unaryMinus
import lib.plus
import lib.minus
import lib.vector.SparseVector
import org.apache.commons.math3.fraction.Fraction

class ABMCMC<AGENT : Agent<AGENT>, ACT: Ordered<ACT>>(
    val model: ABM<AGENT, ACT>,
    nTimesteps: Int,
    val observations: List<Observation<AGENT,ACT>>
) {
    val simplex: SimplexMCMC<Fraction> = SimplexMCMC(
        FractionOperators,
        validTrajectoryConstraints(model, nTimesteps) +
                observations.flatMap { it.eventConstraints() },
        this::logProb
    )


    fun logProb(X: SparseVector<Fraction>): Double {
        val trajectory = Trajectory(model,X)
        val prior = trajectory.logPrior()
        val likelihood = observations.sumByDouble { it.logLikelihood(trajectory) }
        val logP = prior + likelihood
        println("got logprob $prior + $likelihood = $logP")
        return logP
    }


    fun nextSample(): Trajectory<AGENT, ACT> {
        return Trajectory(model, simplex.nextSample())
    }


    companion object {
        fun <AGENT : Agent<AGENT>, ACT: Ordered<ACT>> validTrajectoryConstraints(
            model: ABM<AGENT, ACT>,
            nTimesteps: Int
        ): List<Constraint<Fraction>> {
            println("Constructing model constraints...")
            val constraints = ArrayList<Constraint<Fraction>>(model.actDomain.size*model.agentDomain.size*2)
            constraints.addAll(continuityConstraints(nTimesteps, model))
            constraints.addAll(fermionicConstraints(nTimesteps, model))
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

        fun<AGENT : Agent<AGENT>,ACT: Ordered<ACT>> continuityConstraints(nTimesteps: Int, abm: ABM<AGENT,ACT>): List<Constraint<Fraction>> {
            val nStates = abm.agentDomain.size
//        val acts = abm.actDomain
            val nActs = abm.actDomain.size
            val constraints = ArrayList<Constraint<Fraction>>((nTimesteps-1)*nStates)

            // first do leaving edges
            var i = 0
            for(t in 1 until nTimesteps) {
                for(state in 0 until nStates) {
                    val coeffs = HashMap<Int,Fraction>()
                    val coeffBase = (t*nStates + state)*nActs
                    for(act in 0 until nActs) {
                        coeffs[coeffBase + act] = Fraction(-1)
                    }
                    constraints.add(Constraint(coeffs,"==", Fraction.ZERO))
                }
            }
            // now do incoming edges
            for(state in 0 until nStates) {
                for(act in 0 until nActs) {
                    val consequences = abm.action(abm.agentDomain[state], abm.actDomain[act])
                    for((resultState, n) in consequences) {
                        val resultIndex = resultState.ordinal
                        for(t in 0 until nTimesteps-1) {
                            constraints[t*nStates + resultIndex]
                                .coefficients[(t*nStates + state)*nActs + act] = Fraction(n)
                        }
                    }
                }
            }
//            println("Continuity constraints are $constraints")
            return constraints
        }


        fun<AGENT : Agent<AGENT>,ACT: Ordered<ACT>> fermionicConstraints(nTimesteps: Int, abm: ABM<AGENT,ACT>): List<Constraint<Fraction>> {
            val nStates = abm.agentDomain.size
            val nActs = abm.actDomain.size
            val constraints = ArrayList<Constraint<Fraction>>((nTimesteps+1)*nStates)

            for(t in 0 until nTimesteps) {
                for(state in 0 until nStates) {
                    val coeffs = HashMap<Int,Fraction>()
                    val coeffBase = (t*nStates + state)*nActs
                    for(act in 0 until nActs) {
                        coeffs[coeffBase + act] = Fraction(1)
                    }
                    constraints.add(Constraint(coeffs,"<=", Fraction.ONE))
                }
            }

            val coeffsByState = Array(nStates) { HashMap<Int,Fraction>() }
            for(state in 0 until nStates) {
                val coeffBase = ((nTimesteps-1)*nStates + state)*nActs
                for(act in 0 until nActs) {
                    val consequences = abm.action(abm.agentDomain[state], abm.actDomain[act])
                    for((resultState, n) in consequences) {
                        coeffsByState[resultState.ordinal][coeffBase + act] = Fraction(n)
                    }
                }
            }
            for(state in 0 until nStates) {
                constraints.add(Constraint(coeffsByState[state],"<=", Fraction.ONE))
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
            return if(y.relation == "==") {
                fermionicXImpliesY(x, Constraint(y.coefficients,"<=",y.constant)) +
                        fermionicXImpliesY(x, Constraint(y.coefficients,">=",y.constant))

            } else {
                val const: Fraction
                val coeffs = if (y.relation == "<=") {
                    const = y.constant
                    HashMap(y.coefficients)
                } else {
                    const = -y.constant
                    val negatedCoeffs = HashMap<Int,Fraction>(y.coefficients.size)
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
        fun<AGENT: Agent<AGENT>,ACT: Ordered<ACT>> stateConstraintToActConstraint(stateConstraint: Constraint<Fraction>, timestep: Int, abm: ABM<AGENT,ACT>): Constraint<Fraction> {
            val actCoeffs = HashMap<Int,Fraction>()
            val nActs = abm.actDomain.size
            val timestepBase = abm.agentDomain.size*nActs*timestep
            stateConstraint.coefficients.forEach { (state, coefficient) ->
                for(act in 0 until nActs) {
                    actCoeffs[timestepBase + state*nActs + act] = coefficient
                }
            }
            return Constraint(actCoeffs, stateConstraint.relation, stateConstraint.constant)
        }




    }
}