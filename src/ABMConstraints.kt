import lib.abstractAlgebra.IntOperators
import lib.minus
import lib.plus
import lib.sparseMatrix.GridMapMatrix
import lib.unaryMinus
import org.apache.commons.math3.fraction.Fraction
import java.lang.IllegalArgumentException

object ABMConstraints {
    // nStates      - number of agent states
    // nTimesteps   - number of timesteps
    // nActs        - number of agent actions
    // action       - The "action function", takes an agent state, psi, and an action, a, and returns a sparse vector X
    //                such that x_i is the number of agents in state i after an agent in state psi performs act a.
    fun<AGENT,ACT> continuityConstraints(nTimesteps: Int, abm: ABM<AGENT,ACT>): List<Constraint<Fraction>> {
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
                val consequences = abm.action(abm.agentDomain.toObject(state), abm.actDomain.toObject(act))
                for((resultState, n) in consequences) {
                    val resultIndex = abm.agentDomain.toIndex(resultState)
                    for(t in 0 until nTimesteps-1) {
                        constraints[t*nStates + resultIndex]
                            .coefficients[(t*nStates + state)*nActs + act] = Fraction(n)
                    }
                }
            }
        }
        return constraints
    }


    fun<AGENT,ACT> fermionicConstraints(nTimesteps: Int, abm: ABM<AGENT,ACT>): List<Constraint<Fraction>> {
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
                val consequences = abm.action(abm.agentDomain.toObject(state), abm.actDomain.toObject(act))
                for((resultState, n) in consequences) {
                    coeffsByState[abm.agentDomain.toIndex(resultState)][coeffBase + act] = Fraction(n)
                }
            }
        }
        for(state in 0 until nStates) {
            constraints.add(Constraint(coeffsByState[state],"<=", Fraction.ONE))
        }

        return constraints
    }

    // Returns constraint x -> y
    // under the assumption that
    // 0 <= x <= 1
    // and 0 <= y_i <= 1
    fun fermionicXImpliesY(x: Int, y: Constraint<Fraction>): Constraint<Fraction> {
        if(y.relation == "==") TODO("Not implemented implication to equality, try splitting into two inequalities")
        val coeffs = HashMap(y.coefficients)
        if(y.relation == ">=") coeffs.entries.forEach { it.setValue(-it.value) }
        var maxVal = Fraction.ZERO
        coeffs.values.forEach { if(it > Fraction.ZERO) maxVal += it }
        coeffs[x] = maxVal - y.constant
        return Constraint(coeffs,"<=",maxVal)
    }


    // converts a constraint in terms of state occupation numbers into a constraint on acts
    // in a given timestep
    fun<AGENT,ACT> stateConstraintToActConstraint(stateConstraint: Constraint<Fraction>, timestep: Int, abm: ABM<AGENT,ACT>): Constraint<Fraction> {
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