import lib.*
import java.io.Serializable
import kotlin.random.Random

class PredPreyModel: Hamiltonian<Agent> {
    val deathEvents: Map<Agent,Event<Agent>>

    data class ObservedState(val realState: Multiset<Agent>, val observation: Multiset<Agent>): Serializable

    constructor(params: Params) : super() {
        deathEvents = HashMap()
        for(i in 0 until params.GRIDSIZESQ) {
            Predator(i).hamiltonian(this, params)
            Prey(i).hamiltonian(this, params)
        }
        this.filter { it.consequences.isEmpty() && it.secondaryRequirements.isEmpty() }.forEach {
            deathEvents[it.primaryAgent] = it
        }
    }


    fun sampleContinuousTimePath(startState: Multiset<Agent>, time: Double): List<Pair<Double,Multiset<Agent>>> {
        val distribution = MutableCategorical<Event<Agent>>()
        // initialise distribution
        startState.forEach { distribution.updateRates(it, startState) }

        // initialise path
        val path = ArrayList<Pair<Double,Multiset<Agent>>>((time*distribution.sum()).toInt())
        path.add(Pair(0.0,startState))

        // generate events
        var t = Random.nextExponential(distribution.sum())
        var state = startState
        while(t < time) {
            val nextAct = distribution.sample()
            state = nextAct.actOn(state)
            path.add(Pair(t,state))
            nextAct.modifiedAgents().forEach { distribution.updateRates(it, state) }
            t += Random.nextExponential(distribution.sum())
        }
        return path
    }


    // timestepping
    fun sampleTimesteppingPath(startState: Multiset<Agent>, nSteps: Int): EventTrajectory<Agent> {
        val path = EventTrajectory<Agent>(nSteps)

        // generate events
        var state = startState
        for(t in 1..nSteps) {
            val lastState = state
            val events = ModelEvent<Agent>()
            lastState.forEach { agent ->
                val cumulativeProb = Random.nextDouble()
                val options = primaryRequirementIndex[agent]?:throw(IllegalStateException("No choices for agent $agent"))
                var sum = 0.0
                val nextAct = options.find { act ->
                    sum += act.rateFor(lastState)
                    sum > cumulativeProb
                }?:options.last()
                state = nextAct.actOn(state)
                events.add(nextAct)
            }
            path.add(events)
        }
        return path
    }



    fun generateObservations(startState: Multiset<Agent>, nSteps: Int, pObserve: Double): List<ObservedState> {
        val realTrajectory = sampleTimesteppingPath(startState, nSteps)
        val realStates = realTrajectory.toStateTrajectory()
        val observations = realTrajectory.generateObservations(pObserve)
        return(realStates.zip(observations).map { ObservedState(it.first, it.second) })
    }


    fun MutableCategorical<Event<Agent>>.updateRates(agent: Agent, state: Multiset<Agent>) {
        requirementIndex[agent]?.forEach { act ->
            val rate = act.rateFor(state)
            if(rate == 0.0)
                this.remove(act)
            else
                this[act] = rate
        }
    }

    companion object {
        fun plotTrajectory(trajectory: List<Multiset<Agent>>, gridSize: Int) {
            val data = trajectory.map { frame ->
                frame.supportSet.asSequence().flatMap { agent ->
                    sequenceOf(agent.pos.rem(gridSize), agent.pos.div(gridSize), if(agent is Prey) 1 else 2)
                }.toList()
            }

            gnuplot {
                invoke("set linetype 1 lc 'red'")
                invoke("set linetype 2 lc 'blue'")
                data.forEach { frame ->
                    val pointData = heredoc(frame,3)
                    invoke("plot [0:$gridSize][0:$gridSize] $pointData with points pointtype 5 pointsize 0.5 lc variable")
                    invoke("pause 0.04")
                }
            }
        }


        fun randomState(nAgents: Int, params: Params): Multiset<Agent> {
            val state = HashMultiset<Agent>()
            while(state.size < nAgents) {
                val pos = Random.nextInt(params.GRIDSIZESQ)
                val randAgent = if(Random.nextBoolean()) Predator(pos) else Prey(pos)
                state.add(randAgent)
            }
            return state
        }


        fun randomState(nPredator: Int, nPrey: Int, params: Params): Multiset<Agent> {
            val state = HashMultiset<Agent>()
            for(i in 1..nPredator) {
                val pos = Random.nextInt(params.GRIDSIZESQ)
                state.add(Predator(pos))
            }
            for(i in 1..nPrey) {
                val pos = Random.nextInt(params.GRIDSIZESQ)
                state.add(Prey(pos))
            }
            return state
        }

    }
}