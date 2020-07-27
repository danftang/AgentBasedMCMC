import com.google.ortools.linearsolver.MPSolver
import lib.Multiset
import lib.emptyMultiset
import java.io.*
import java.lang.IllegalStateException
import java.util.*
import kotlin.system.measureTimeMillis
import kotlin.time.measureTime

class MAPOrbitSolver<AGENT>: Serializable {
    val timesteps = ArrayList<Timestep<AGENT>>()
    val startState: StartState<AGENT>
    val hamiltonian: Hamiltonian<AGENT>
//    var solver: MPSolver

    val trajectory: EventTrajectory<AGENT>
        get() {
            val result = EventTrajectory<AGENT>(timesteps.size)
            timesteps.mapTo(result) { it.committedEvents }
            return result
        }


    val observations: List<Multiset<AGENT>>
        get() = timesteps.map { it.observations }



    constructor(hamiltonian: Hamiltonian<AGENT>, startState: Multiset<AGENT>) {
        this.hamiltonian = hamiltonian
        this.startState = StartState(startState)
    }


    fun addObservation(observation: Multiset<AGENT>) {
        val lastState: UnknownModelState<AGENT> = if(timesteps.isEmpty()) startState else timesteps.last()
        timesteps.add(Timestep(hamiltonian, lastState, observation))
    }


    fun solutionIsCorrect(isPartial: Boolean): Boolean {
        if(!trajectory.isFeasible(startState.committedConsequences, isPartial)) return false
        timesteps.forEach { step ->
            if(!step.observations.isSubsetOf(step.committedEvents.consequences)) return false
        }
        return true
    }


    fun addObservations(observations: Iterable<Multiset<AGENT>>) {
        observations.forEach { addObservation(it) }
    }


    fun completeSolve() {
        println("Starting complete solve")
        val solver = MPSolver("HamiltonianSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        timesteps.forEach { it.addForwardEvents() }

        timesteps.forEach { it.setupProblem(solver, true) }
        println("solving for ${solver.numVariables()} variables and ${solver.numConstraints()} constraints")
            val solveState = solver.solve()
            println("solveState = $solveState")

        timesteps.forEach {
            it.applySolution()
        }
        startState.endStateBinaryIndicators.clear()
        startState.endStateIndicators.clear()

    }

    fun minimalSolve() {
        do {
            val solver = MPSolver("HamiltonianSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
            calculateForwardBackwardPotentialEvents()

            timesteps.forEach { it.setupProblem(solver, false) }
            println("solving for ${solver.numVariables()} variables and ${solver.numConstraints()} constraints")
            val solveState = solver.solve()
            println("solveState = $solveState")

            if (solveState == MPSolver.ResultStatus.INFEASIBLE) {
                println("COULDN'T SOLVE. ROLLING BACK...")
                val timestepIt = timesteps.asReversed().iterator()
                var doneRollback = false
                while(timestepIt.hasNext() && !doneRollback) {
                    doneRollback = timestepIt.next().rollback()
                }
                if(!doneRollback) throw(IllegalStateException("Unsolvable state"))
                timesteps.forEach { it.clearAllIndicators() }
            } else {
                timesteps.forEach { it.applySolution() }
            }
            startState.endStateBinaryIndicators.clear()
            startState.endStateIndicators.clear()

        } while(solveState == MPSolver.ResultStatus.INFEASIBLE)

    }


    fun removeDeadAgents(maxStepsUnseen: Int) {
        timesteps.asReversed().asSequence().drop(maxStepsUnseen-1).forEach { timestep ->
            timestep.hangingAgents.forEach { agent ->
                timestep.committedEvents.add(Event(setOf(agent), emptySet(), emptyMultiset(), 1.0, agent))
            }
        }
    }



    private fun calculateForwardBackwardPotentialEvents() {
        timesteps.forEach { it.addForwardEvents() }
        var forwardRequirements = emptySet<AGENT>()
        var forwardCommittedAbsences = emptySet<AGENT>()
        timesteps.asReversed().forEach { timestep ->
            timestep.filterBackwardEvents(forwardRequirements, forwardCommittedAbsences)
            forwardRequirements = timestep.potentialRequirementsFootprint
            forwardCommittedAbsences = timestep.committedAbsenceFootprint
        }
    }


}