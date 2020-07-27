import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import lib.HashMultiset
import lib.Multiset
import lib.MutableMultiset
import java.io.Serializable

class StartState<AGENT>: UnknownModelState<AGENT>, Serializable {
    override val committedConsequences: MutableMultiset<AGENT>
    override val potentialConsequencesFootprint: Set<AGENT>
        get() = emptySet()

    val endStateIndicators = HashMap<AGENT,MPVariable>()
    val endStateBinaryIndicators = HashMap<AGENT,MPVariable>()

    constructor(state: Multiset<AGENT>) {
        committedConsequences = HashMultiset()
        committedConsequences.addAll(state)
    }

    override fun getEndStateVriable(solver: MPSolver, agent: AGENT): MPVariable {
        return endStateIndicators.getOrPut(agent) {
            val multiplicity = committedConsequences.count(agent).toDouble()
            solver.makeIntVar(multiplicity, multiplicity,  Timestep.getNextVarId())
        }
    }

    override fun getEndStateBinaryVriable(solver: MPSolver, agent: AGENT): MPVariable {
        return endStateBinaryIndicators.getOrPut(agent) {
            val isPresent = if(committedConsequences.contains(agent)) 1.0 else 0.0
            solver.makeIntVar(isPresent, isPresent,  Timestep.getNextVarId())
        }
    }

//    override fun addConsequencesToConstraints(constraints: Map<AGENT, MPConstraint>, multiplier: Double) {
//    }

//    override fun hasIndicators() = false

}