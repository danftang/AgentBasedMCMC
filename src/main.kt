fun main() {
    System.loadLibrary("jniortools")

    val WINDOW_LEN = 5
    val TOTAL_STEPS = 15
    val params = StandardParams
    val myModel = PredPreyModel(params)
    val startState = PredPreyModel.randomState(50, params)
    val observations = myModel.generateObservations(startState, TOTAL_STEPS, 1.0)
    val windows = observations
        .drop(1)
        .windowed(WINDOW_LEN, WINDOW_LEN, true) {window ->
            window.map { obs -> obs.observation }
        }
    println("Real orbit")
    observations.forEach {println(it.realState)}
    println("Observations")
    observations.forEach {println(it.observation)}



    val mySolver = MAPOrbitSolver(myModel, startState)

    windows.forEach {window ->
        println("Adding window $window")
        mySolver.addObservations(window)
        mySolver.minimalSolve()

        println("MAP orbit is")
        mySolver.timesteps.forEach { println(it.committedEvents) }
        println("history is")
        println(startState)
        mySolver.timesteps.forEach { println(it.committedConsequences) }

//        removeDeadAgents(mySolver, myModel, 8)
    }
}

//fun removeDeadAgents(solver: MAPOrbitSolver<Agent>, model: PredPreyModel, maxStepsUnseen: Int) {
//    solver.timesteps.descendingIterator().asSequence().drop(maxStepsUnseen-1).forEach { timestep ->
//        timestep.previousState.sources.forEach { agent ->
//            model.deathEvents[agent]?.also {
//                timestep.committedEvents.add(it)
//            }
//        }
//        timestep.previousState.sources.clear()
//    }
//}
