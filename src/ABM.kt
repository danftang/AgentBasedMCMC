import lib.collections.Multiset
import lib.gnuplot
import org.apache.commons.math3.fraction.Fraction
import kotlin.math.ln
import kotlin.random.Random

interface ABM<AGENT: Agent<AGENT>,ACT: Ordered<ACT>> {

    data class Event<AGENT: Ordered<AGENT>,ACT: Ordered<ACT>>(val time: Int, val agent: AGENT, val act: ACT): Ordered<Event<AGENT,ACT>> {
        override val ordinal: Int
            get() = time * agent.domain.size * act.domain.size +
                    agent.ordinal*act.domain.size +
                    act.ordinal

        override val domain: CountableDomain<Event<AGENT, ACT>>
            get() = domain(agent.domain,act.domain)

        companion object {
            fun <AGENT : Ordered<AGENT>, ACT : Ordered<ACT>> domain(
                agentDomain: CountableDomain<AGENT>,
                actDomain: CountableDomain<ACT>
            ): CountableDomain<Event<AGENT, ACT>> =
                object : CountableDomain<Event<AGENT, ACT>> {
                    override val size: Int
                        get() = Int.MAX_VALUE

                    override fun get(index: Int): Event<AGENT, ACT> {
                        return Event(
                            index / (agentDomain.size * actDomain.size),
                            agentDomain.get(index.rem(agentDomain.size * actDomain.size) / actDomain.size),
                            actDomain.get(index.rem(actDomain.size))
                        )
                    }
                }
        }
    }

    val actDomain: CountableDomain<ACT>
    val agentDomain: CountableDomain<AGENT>
    val eventDomain: CountableDomain<Event<AGENT,ACT>>
        get() = Event.domain(agentDomain, actDomain)

    fun isFermionic(): Boolean = true
    fun action(startState: AGENT, act: ACT): Map<AGENT,Int>
    fun timestepSupport(agent: AGENT, act: ACT): List<Constraint<Fraction>> // to be generated automatically...eventually.


    fun runABM(startState: Multiset<AGENT>, nTimesteps: Int): Trajectory<AGENT, ACT> {
        var state = startState
        val trajectory = Trajectory<AGENT,ACT>()
        for(t in 0 until nTimesteps) {
            for((agent, occupation) in state.nonZeroEntries) {
                for(i in 1..occupation) {
                    val act = agent.concreteTimestep(state)
                    trajectory[t,agent,act] += 1
                }
            }
        }
        return trajectory
    }


    // Assumes trajectory is valid
    fun logProb(trajectory: Trajectory<AGENT, ACT>): Double {
        var logProb = 0.0
        for(time in trajectory.indices) {
            val state = trajectory.stateAt(time)
            for((agentAct, occupation) in trajectory[time].nonZeroEntries) {
                val (agent,act) = agentAct
                logProb += occupation*ln(agent.timestep(state)[act.ordinal])
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
        return actDomain.get(actIndex)
    }

    fun Multiset<Pair<AGENT,ACT>>.toABMState(): Multiset<AGENT> {
        val state = Multiset<AGENT>()
        for((agentAct, occupation) in nonZeroEntries) {
            val (agent, _) = agentAct
            state[agent] += occupation
        }
        return state
    }

    fun Trajectory<AGENT, ACT>.nAgents(time: Int, agent: AGENT): Int {
        var count = 0
        for(act in 0 until actDomain.size) {
            if(this[time][Pair(agent,actDomain.get(act))] >= 1) count++
        }
        return count
    }

    // Plots a Feynmann diagram of this trajectory, time on the y axis
    // agent state on the x axis and vectors for acts
    fun plot(trajectory: Trajectory<AGENT, ACT>) {
        var maxX = 0
        var maxY = 0
        val particleVectors = ArrayList<Int>()      // in format (x,y,dx,dy,colour)...
        for ((event, occupation) in trajectory.events) {
            val consequences = action(event.agent, event.act)
            val x0 = event.agent.ordinal
            val y0 = event.time
            if(y0 >= maxY) maxY = y0+1
            if(x0 > maxX) maxX = x0
            for ((endPoint, multiplicity) in consequences) {
                val x1 = endPoint.ordinal
                if(x1 > maxX) maxX = x1
                particleVectors.addAll(arrayOf(x0, y0, x1 - x0, 1, if(occupation < 0) 0xff0000 else 0xff ))
            }
        }

        gnuplot {
            if (particleVectors.size > 0) {
                val particleData = heredoc(particleVectors, 5)
                invoke("plot [0:${maxX}][0:$maxY] $particleData with vectors lw 2 lc rgbcolor variable")
            }
        }
    }

}