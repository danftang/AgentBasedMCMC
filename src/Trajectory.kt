import lib.collections.Multiset
import lib.vector.SparseVector
import org.apache.commons.math3.fraction.Fraction
import java.lang.StringBuilder
import java.util.AbstractMap
import kotlin.math.ln

class Trajectory<AGENT : Agent<AGENT>, ACT: Ordered<ACT>>(val timesteps: ArrayList<Multiset<Pair<AGENT, ACT>>> = ArrayList()) :
    MutableList<Multiset<Pair<AGENT, ACT>>> by timesteps {

    val events = timesteps.asSequence()
        .mapIndexed { time, step ->
            step.nonZeroEntries.asSequence().map { (agentAct, occupation) ->
                AbstractMap.SimpleEntry(ABM.Event(time, agentAct.first, agentAct.second), occupation)
            }
        }
        .flatten()

    constructor(model: ABM<AGENT,ACT>, actVector: SparseVector<Fraction>): this() {
        for((actId, occupation) in actVector.nonZeroEntries) {
            this[model.eventDomain[actId]] = occupation.toInt()
        }
    }

    inline operator fun get(time: Int, agent: AGENT, act: ACT): Int {
        return if (time < timesteps.size) timesteps[time][Pair(agent, act)] else 0
    }

    inline operator fun get(event: ABM.Event<AGENT, ACT>): Int {
        return get(event.time, event.agent, event.act)
    }

    inline operator fun set(time: Int, agent: AGENT, act: ACT, occupation: Int) {
        while (time >= timesteps.size) timesteps.add(Multiset())
        timesteps[time][Pair(agent, act)] = occupation
    }

    inline operator fun set(event: ABM.Event<AGENT, ACT>, occupation: Int) {
        set(event.time, event.agent, event.act, occupation)
    }

    fun stateAt(time: Int): Multiset<AGENT> {
        val state = Multiset<AGENT>()
        for((agentAct, occupation) in this[time].nonZeroEntries) {
            val (agent, _) = agentAct
            state[agent] += occupation
        }
        return state
    }

    fun logPrior(): Double {
        var logP = 0.0
        for(t in indices) {
            val state = stateAt(t)
            for((agentAct, occupation) in this[t].nonZeroEntries) {
                val (agent, act) = agentAct
                logP += occupation*ln(agent.timestep(state)[act.ordinal])
            }
        }
        return logP
    }


    override fun toString(): String {
        val out = StringBuilder()
        for(step in timesteps) {
            out.append(step.toString())
            out.append('\n')
        }
        return out.toString()
    }


}
