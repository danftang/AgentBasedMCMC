import lib.collections.Multiset
import lib.vector.SparseVector
import org.apache.commons.math3.fraction.Fraction
import kotlin.math.ln
import kotlin.random.Random

interface ABM<AGENT: Agent<AGENT>,ACT> {

    class Trajectory<AGENT:Agent<AGENT>,ACT>(val timesteps: ArrayList<Multiset<Pair<AGENT,ACT>>> = ArrayList()):
        MutableList<Multiset<Pair<AGENT,ACT>>> by timesteps {
        }

    val actDomain: CountableDomain<ACT>
    val agentDomain: CountableDomain<AGENT>

    fun isFermionic(): Boolean = true
    fun trajectoryLogProb(x: SparseVector<Fraction>): Double = 0.0
    fun action(startState: AGENT, act: ACT): Map<AGENT,Int>
    fun timestepSupport(agent: AGENT, act: ACT): List<Constraint<Fraction>> // to be generated automatically...eventually.


    fun runABM(startState: Multiset<AGENT>, nTimesteps: Int): Trajectory<AGENT,ACT> {
        var state = startState
        val trajectory = Trajectory<AGENT,ACT>()
        for(t in 0 until nTimesteps) {
            val timestep = Multiset<Pair<AGENT,ACT>>()
            for((agent, occupation) in state) {
                for(i in 1..occupation) {
                    val act = agent.concreteTimestep(state)
                    timestep[Pair(agent,act)] += 1
                }
            }
            trajectory.add(timestep)
        }
        return trajectory
    }


    // Assumes trajectory is valid
    fun logProb(trajectory: Trajectory<AGENT,ACT>): Double {
        var logProb = 0.0
        for(timestep in trajectory) {
            val state = timestep.toABMState()
            for((agentAct, occupation) in timestep) {
                val (agent,act) = agentAct
                logProb += occupation*ln(agent.timestep(state)[actDomain.toIndex(act)])
            }
        }
        return logProb
    }


    fun AGENT.concreteTimestep(others: Multiset<AGENT>): ACT {
        val actDistribution = this.timestep(others)
        val r = Random.nextDouble()
        var actIndex = 0
        var cumulativeDistribution = actDistribution[0]
        while(r > cumulativeDistribution) {
            cumulativeDistribution += actDistribution[++actIndex]
        }
        return actDomain.toObject(actIndex)
    }

    fun Multiset<Pair<AGENT,ACT>>.toABMState(): Multiset<AGENT> {
        val state = Multiset<AGENT>()
        for((agentAct, occupation) in this) {
            val (agent, _) = agentAct
            state[agent] += occupation
        }
        return state
    }

    fun Trajectory<AGENT,ACT>.nAgents(agent: AGENT, time: Int): Int {
        var count = 0
        for(act in 0 until actDomain.size) {
            if(this[time][Pair(agent,actDomain.toObject(act))] >= 1) count++
        }
        return count
    }
}