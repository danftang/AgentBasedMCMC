import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import com.google.ortools.sat.CpModel
import com.google.ortools.sat.IntVar
import lib.*
import java.io.Serializable
import kotlin.math.ln

class Timestep<AGENT>: UnknownModelState<AGENT>, Serializable {
    val hamiltonian: Hamiltonian<AGENT>
    val previousState: UnknownModelState<AGENT>
    val observations = HashMultiset<AGENT>()
    var potentialEvents: Set<Event<AGENT>> = setOf()
    val committedEvents = ModelEvent<AGENT>()
    val eventIndicators = HashMap<Event<AGENT>, IntVar>()
    val endStateIndicators = HashMap<AGENT, IntVar>()
    val endStateBinaryIndicators = HashMap<AGENT,IntVar>()
    val M = 100.0 // Maximum multiplicity


    override val committedConsequences: Multiset<AGENT>
        get() = committedEvents.consequences()

    override val potentialConsequencesFootprint: Set<AGENT>
        get() = potentialEvents.asSequence().flatMap { it.consequences.supportSet.asSequence() }.toSet()

    val potentialRequirementsFootprint: Set<AGENT>
        get() = potentialEvents.asSequence().flatMap { it.presenceRequirements.asSequence() }.toSet()

    val committedAbsenceFootprint: Set<AGENT>
        get() = committedEvents.asSequence().flatMap { it.absenceRequirements.asSequence() }.toSet()

    val hangingAgents: Multiset<AGENT>
        get() = previousState.committedConsequences - committedEvents.primaryRequirements()

//    val absenceFootprint: Set<AGENT>
//        get() = potentialEvents.asSequence().flatMap { it.absenceRequirements.asSequence() }.toSet()



    constructor(hamiltonian: Hamiltonian<AGENT>, previousState: UnknownModelState<AGENT>, observations: Multiset<AGENT>) {
        this.hamiltonian = hamiltonian
        this.previousState = previousState
        this.observations.addAll(observations)
    }



    fun addForwardEvents() {
        val previousPotentialConsequences = previousState.potentialConsequencesFootprint
        val previousCommittedFootprint = previousState.committedConsequences.supportSet
        val potentialPrimaryAgents = previousPotentialConsequences + hangingAgents.supportSet
        val potentialPresentAgents = previousPotentialConsequences + previousCommittedFootprint

        potentialEvents = hamiltonian
            .eventsPresenceSatisfiedBy(potentialPresentAgents)
            .asSequence()
            .filter { potentialPrimaryAgents.contains(it.primaryAgent) }
            .filter { it.absenceRequirements.isDisjoint(previousCommittedFootprint) }
            .toSet()
    }



    fun observationsAreSatisfied(): Boolean = committedEvents.isNotEmpty()


    fun filterBackwardEvents(nextRequirementsFootprint: Set<AGENT>, nextCommittedAbsences: Set<AGENT>) {
        val allRequirements = if(observationsAreSatisfied())
            nextRequirementsFootprint
        else
            nextRequirementsFootprint.union(observations.supportSet)
        potentialEvents = potentialEvents.asSequence().filter { event ->
            event.consequences.supportSet.isDisjoint(nextCommittedAbsences) &&
                    !event.consequences.supportSet.isDisjoint(allRequirements)
        }.toSet()
    }


    fun setupProblem(model: CpModel, isComplete: Boolean) {
        // setup vars
        eventIndicators.clear()
        endStateIndicators.clear()
        endStateBinaryIndicators.clear()
        potentialEvents.forEach {event ->
            eventIndicators[event] = model.newIntVar(0, Long.MAX_VALUE, getNextVarId())
        }

        setupObservationConstraints(model)

        val primaryConstraints = setupRequirementsConstraints(model, -1.0) { setOf(it.primaryAgent) }
        val primaryCommittments = committedEvents.primaryRequirements()
        primaryConstraints.forEach { (agent, constraint) ->
            constraint.setCoefficient(previousState.getEndStateVriable(model, agent), 1.0)
            val primaryCount = primaryCommittments.count(agent).toDouble()
            if(isComplete)
                constraint.setBounds(primaryCount, primaryCount)
            else
                constraint.setBounds(primaryCount,Double.POSITIVE_INFINITY)
        }

        val secondaryConstraints = setupRequirementsConstraints(model, -1.0) { it.secondaryRequirements }
        secondaryConstraints.forEach { (agent, constraint) ->
            constraint.setCoefficient(previousState.getEndStateBinaryVriable(model, agent), M)
            constraint.setBounds(0.0, Double.POSITIVE_INFINITY)
        }

        val absenceConstraints = setupRequirementsConstraints(model, 1.0) { it.absenceRequirements }
        absenceConstraints.forEach { (agent, constraint) ->
            constraint.setCoefficient(previousState.getEndStateBinaryVriable(model, agent), M)
            constraint.setBounds(0.0, if(committedAbsenceFootprint.contains(agent)) M-1.0 else M)
        }

        setupObjective(model)
    }


    fun clearAllIndicators() {
        eventIndicators.clear()
        endStateBinaryIndicators.clear()
        endStateIndicators.clear()
    }


    fun applySolution() {
        val newCommitedEvents = ArrayList<Event<AGENT>>()
        eventIndicators.forEach { (event, mpVar) ->
            for(i in 1..mpVar.solutionValue().toInt()) newCommitedEvents.add(event)
        }
        clearAllIndicators()
        committedEvents.addAll(newCommitedEvents)

        // check secondary requirements are met
        newCommitedEvents.forEach {
            if(!previousState.committedConsequences.containsAll(it.secondaryRequirements))
                throw(IllegalStateException("Solution orbit does not meet secondary requirements for event $it"))
        }

        // check absence requirements are met
        newCommitedEvents.forEach {
            if(!previousState.committedConsequences.intersect(it.absenceRequirements).isEmpty())
                throw(IllegalStateException("Solution orbit does not meet absence requirements for event $it"))
        }

        // check primary requirements are met
        if(!committedEvents.primaryRequirements().isSubsetOf(previousState.committedConsequences))
            throw(IllegalStateException("Solution orbit does not meet primary requirements"))

        // check observations are met
        if(!observations.isSubsetOf(committedConsequences)) throw(IllegalStateException("Solution orbit does not satisfy observations."))
    }



    fun setupObservationConstraints(solver: CpModel) {
        if(observationsAreSatisfied()) return
        for ((agent, nObserved) in observations.counts) {
            val satisfyObservationConstraint = solver.makeConstraint(nObserved.toDouble(), Double.POSITIVE_INFINITY)
            hamiltonian.eventsWithConsequence(agent).forEach { event ->
                eventIndicators[event]?.also { indicatorVar ->
                    satisfyObservationConstraint.setCoefficient(indicatorVar, event.consequences.count(agent).toDouble())
                }
            }
        }
    }


    // Creates a set of constraints such that for each member of requiredAgents for each potential event,
    // the following coefficients are added to the constraint
    // multiplier * (sum_{event_t} agent \in requiredAgents(event)*i_event)
    fun setupRequirementsConstraints(solver: CpModel, multiplier: Double, requiredAgents: (Event<AGENT>) -> Set<AGENT>): Map<AGENT, MPConstraint> {
        val requirementConstraints = HashMap<AGENT, MPConstraint>()
        eventIndicators.forEach { (event, mpvar) ->
            requiredAgents(event).forEach { agent ->
                val constraint = requirementConstraints.getOrPut(agent) { solver.makeConstraint() }
                constraint.setCoefficient(mpvar, multiplier)
            }
        }
        return requirementConstraints
    }


    fun setupObjective(solver: MPSolver) {
        val obj = solver.objective()
        eventIndicators.forEach { (event, mpVar) ->
            obj.setCoefficient(mpVar, ln(event.rate))
        }
        obj.setMaximization()
    }




    fun rollback(): Boolean {
        if(committedEvents.isEmpty()) return false
        committedEvents.clear()
        return true
    }


    companion object {
        var mpVarId = 0

        fun getNextVarId(): String = mpVarId++.toString()
    }

    override fun getEndStateVriable(solver: CpModel, agent: AGENT): IntVar {
        return endStateIndicators.getOrPut(agent) {
            val commitment = committedConsequences.count(agent).toDouble()
            val indicator = solver.makeIntVar(0.0, Double.POSITIVE_INFINITY,  getNextVarId())
            val constraint = solver.makeConstraint(commitment, commitment)
            constraint.setCoefficient(indicator, 1.0)
            hamiltonian.eventsWithConsequence(agent).forEach { event ->
                eventIndicators[event]?.also { eventIndicator ->
                    constraint.setCoefficient(eventIndicator, -event.consequences.count(agent).toDouble())
                }
            }
            indicator
        }
    }

    override fun getEndStateBinaryVriable(solver: CpModel, agent: AGENT): IntVar {
        return endStateBinaryIndicators.getOrPut(agent) {
            val psi_agent = getEndStateVriable(solver, agent)
            val binaryIndicator = solver.makeIntVar(0.0, 1.0, getNextVarId())
            val gtConstraint = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY)
            gtConstraint.setCoefficient(binaryIndicator, M)
            gtConstraint.setCoefficient(psi_agent, -1.0)
            val ltConstraint = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY)
            ltConstraint.setCoefficient(binaryIndicator, -1.0)
            ltConstraint.setCoefficient(psi_agent, 1.0)
            binaryIndicator
        }
    }


}