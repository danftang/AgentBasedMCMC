import lib.Multiset
import lib.divergence
import java.io.FileOutputStream
import java.io.ObjectOutputStream
import java.io.Serializable
import kotlin.system.measureTimeMillis

class PredPreyProblem : Serializable {

    val realTrajectory: EventTrajectory<Agent>
    val solver: MAPOrbitSolver<Agent>
    val params: Params

    constructor(realTrajectory: EventTrajectory<Agent>, solver: MAPOrbitSolver<Agent>, params: Params) {
        this.realTrajectory = realTrajectory
        this.solver = solver
        this.params = params
    }

    constructor(realTrajectory: EventTrajectory<Agent>, params: Params): this(
        realTrajectory,
        MAPOrbitSolver(PredPreyModel(params), realTrajectory.impliedStartState()),
        params
    )

    constructor(
        nPredator: Int,
        nPrey: Int,
        nSteps: Int,
        params: Params
    ) {
        this.params = params
        val myModel = PredPreyModel(params)
        val startState = PredPreyModel.randomState(nPredator, nPrey, params)
        realTrajectory = PredPreyModel(params).sampleTimesteppingPath(startState, nSteps)
        solver = MAPOrbitSolver(myModel, startState)
//        solver.addObservations(realTrajectory.generateObservations(pObserve))
    }


    fun onlineSolve(observations: List<Multiset<Agent>>, windowLen: Int, removeDeadAfter: Int) {
        val windows = observations.windowed(windowLen, windowLen, true)

            for (i in windows.indices) {
                val window = windows[i]
                println("Adding window $i of ${windows.size - 1}")
                solver.addObservations(window)
                solver.minimalSolve()
                solver.removeDeadAgents(removeDeadAfter)
            }
            println("Completing trajectory")
            solver.completeSolve()
    }

    fun offlineSolve(observations: List<Multiset<Agent>>) {
        solver.addObservations(observations)
        solver.completeSolve()
    }

    fun calcDivergences(): Pair<List<Double>, List<Double>> {
        val unobservedPosterior = solver.observations
            .zip(solver.trajectory)
            .map {(observation, posterior) ->
                posterior.consequences - observation
            }

        val unobservedReal = realTrajectory.toStateTrajectory()
            .zip(solver.observations)
            .map {(real, obs) ->
                real - obs
            }

        val divergences = unobservedReal
            .zip(unobservedPosterior)
            .map { (real, predicted) ->
                predicted.divergence(real, params.GRIDSIZE)
            }

        val refDivergences = unobservedReal
            .map { realState ->
                PredPreyModel.randomState(100,params).divergence(realState, params.GRIDSIZE)
            }
        return Pair(divergences, refDivergences)
    }

    companion object {
        private const val serialVersionUID: Long = 7222081338134265798
    }
}