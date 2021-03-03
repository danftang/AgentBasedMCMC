import lib.collections.Multiset
import lib.gnuplot
import org.apache.commons.math3.fraction.Fraction
import java.lang.RuntimeException
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
    fun consequences(startState: AGENT, act: ACT): Multiset<AGENT>

    // returns the constraints implied by the given event
    fun timestepEventConstraints(event: Event<AGENT,ACT>): List<Constraint<Fraction>> // to be generated automatically...eventually.


    fun fermionicRunABM(startState: Multiset<AGENT>, nTimesteps: Int): Trajectory<AGENT, ACT> {
        var t0State = startState
        var t1State: Multiset<AGENT>
        val trajectory = Trajectory(this)
        for(t in 0 until nTimesteps) {
            t1State = Multiset()
            for((agent, occupation) in t0State.entries) {
                assert(occupation == 1)
                val act = agent.fermionicConcreteTimestep(t0State, t1State)
                trajectory[t,agent,act] = 1
                t1State += consequences(agent,act)
            }
            t0State = t1State
        }
        return trajectory
    }


    fun nonFermionicRunABM(startState: Multiset<AGENT>, nTimesteps: Int): Trajectory<AGENT, ACT> {
        var t0State = startState
        var t1State: Multiset<AGENT>
        val trajectory = Trajectory(this)
        for(t in 0 until nTimesteps) {
            t1State = Multiset()
            for((agent, occupation) in t0State.entries) {
                for(i in 1..occupation) {
                    val act = agent.nonFermionicConcreteTimestep(t0State)
                    trajectory[t,agent,act] += 1
                    t1State += consequences(agent,act)
                }
            }
            t0State = t1State
        }
        return trajectory
    }


    // Assumes trajectory is valid
    fun logProb(trajectory: Trajectory<AGENT, ACT>): Double {
        var logProb = 0.0
        for(time in 0 until trajectory.nTimesteps) {
            val state = trajectory.stateAt(time)
            for((agentAct, occupation) in trajectory[time].entries) {
                val (agent,act) = agentAct
                logProb += occupation*ln(agent.timestep(state)[act.ordinal])
            }
        }
        return logProb
    }


    // Fermionic timestep
    // this is a bit of a quick and dirty method by
    // just choosing an act that doesn't overlap earlier agents
    // but should work when death is an action
    fun AGENT.fermionicConcreteTimestep(others: Multiset<AGENT>, endStateOccupation: Multiset<AGENT>): ACT {
        val actDistribution = this.timestep(others)
        val fermionicActDistribution = actDomain
            .filter { act ->
                consequences(this, act).entries.all { endStateOccupation[it.key] + it.value <= 1 }
            }
            .map { Pair(it, actDistribution[it.ordinal]) }
        val sumOfProbs = fermionicActDistribution.sumByDouble { it.second }
        val r = Random.nextDouble()*sumOfProbs
        var cumulativeDistribution = 0.0
        return fermionicActDistribution
            .find {
                cumulativeDistribution += it.second
                r < cumulativeDistribution
            }
            ?.first
            ?:throw(RuntimeException("No fermionic acts available"))
    }

    fun AGENT.nonFermionicConcreteTimestep(others: Multiset<AGENT>): ACT {
        val actDistribution = this.timestep(others)
        val r = Random.nextDouble()
        var actIndex = 0
        var cumulativeDistribution = actDistribution[0]
        while(r > cumulativeDistribution) {
            cumulativeDistribution += actDistribution[++actIndex]
        }
        return actDomain.get(actIndex)
    }


    // converts a constraint in terms of state occupation numbers into a constraint on acts
    // in a given timestep
//    fun stateToActConstraint(stateConstraint: Constraint<Fraction>, timestep: Int) = stateConstraint.stateToActConstraint(timestep)
    fun Constraint<Fraction>.stateToEventConstraint(timestep: Int): Constraint<Fraction> {
        val actCoeffs = HashMap<Int,Fraction>()
        val nActs = actDomain.size
        val timestepBase = agentDomain.size*nActs*timestep
        coefficients.forEach { (state, coefficient) ->
            for(act in 0 until nActs) {
                actCoeffs[timestepBase + state*nActs + act] = coefficient
            }
        }
        return Constraint(actCoeffs, relation, constant)
    }



    // Plots a Feynmann diagram of this trajectory, time on the y axis
    // agent state on the x axis and vectors for acts
    fun plot(trajectory: Trajectory<AGENT, ACT>) {
        var maxX = 1
        var maxY = 1
        val particleVectors = ArrayList<Int>()      // in format (x,y,dx,dy,colour)...
        for ((event, occupation) in trajectory.events) {
            val consequences = consequences(event.agent, event.act)
            val x0 = event.agent.ordinal
            val y0 = event.time
            if(y0 >= maxY) maxY = y0+1
            if(x0 > maxX) maxX = x0
            for ((endPoint, multiplicity) in consequences.entries) {
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