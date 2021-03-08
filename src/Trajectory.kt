import lib.abstractAlgebra.FractionOperators
import lib.collections.Multiset
import lib.collections.emptyMultiset
import lib.sparseVector.SparseVector
import lib.sparseVector.asVector
import org.apache.commons.math3.fraction.Fraction
import org.apache.commons.math3.util.CombinatoricsUtils
import java.lang.StringBuilder
import java.util.AbstractMap
import kotlin.math.ln
import kotlin.math.roundToInt

class Trajectory<AGENT : Agent<AGENT>, ACT: Ordered<ACT>>(
    val model: ABM<AGENT, ACT>,
    val eventVector: SparseVector<Fraction>
) {

    val stateTrajectory: List<Multiset<AGENT>> by lazy {
        val stateTrajectory = ArrayList<Multiset<AGENT>>()
        for((event, occupation) in events) {
            while(stateTrajectory.size <= event.time) stateTrajectory.add(Multiset())
            stateTrajectory[event.time][event.agent] += occupation
        }
        stateTrajectory
    }

    val finalState: Multiset<AGENT> by lazy {
        val finalState = Multiset<AGENT>()
        val penultimateTime = stateTrajectory.lastIndex
        stateTrajectory
            .lastOrNull()
            ?.members
            ?.forEach { agent ->
                for(act in model.actDomain) {
                    val occupation = this[penultimateTime, agent, act]
                    if(occupation > 0) {
                        finalState += model.consequences(agent, act) * occupation
                    }
                }
            }
        finalState
    }


    val events: Sequence<AbstractMap.SimpleEntry<ABM.Event<AGENT, ACT>, Int>>
        get() = eventVector.nonZeroEntries
            .asSequence()
            .map { AbstractMap.SimpleEntry(model.eventDomain[it.key], it.value.toDouble().roundToInt()) }


    inline operator fun get(time: Int, agent: AGENT, act: ACT): Int {
        return get(ABM.Event(time, agent, act))
    }

    inline operator fun get(event: ABM.Event<AGENT, ACT>): Int {
        return eventVector[event.ordinal].toDouble().roundToInt()
    }



    fun stateAt(time: Int): Multiset<AGENT> {
        if(time < stateTrajectory.size) return stateTrajectory[time]
        if(time == stateTrajectory.size && time != 0) return finalState
        return emptyMultiset()
    }


    fun nAgents(time: Int, agent: AGENT): Int {
        if(time < stateTrajectory.size) return stateTrajectory[time][agent]
        if(time == stateTrajectory.size && time != 0) return finalState[agent]
        return 0
    }


    fun logPrior(): Double {
        var logP = 0.0
        for((event,occupation) in events) {
            logP += occupation * ln( event.agent.timestep(stateTrajectory[event.time])[event.act.ordinal] )
//                   - CombinatoricsUtils.factorialLog(occupation) // for non-fermionic
        }
// Add this for non-fermionic
//        for(state in stateTrajectory) {
//            for((_, occupation) in state.entries) {
//                logP += CombinatoricsUtils.factorialLog(occupation)
//            }
//        }
        return logP
    }


    override fun toString(): String {
        val out = StringBuilder()
        for(step in events) {
            out.append(step.toString())
            out.append('\n')
        }
        return out.toString()
    }



//    val model: ABM<AGENT,ACT>
//
//    val timesteps: ArrayList<Multiset<Pair<AGENT, ACT>>> = ArrayList()
//
//    val nTimesteps: Int
//        get() = timesteps.size
//
//    val events: Sequence<AbstractMap.SimpleEntry<ABM.Event<AGENT, ACT>, Int>>
//        get() = timesteps.asSequence()
//            .mapIndexed { time, step ->
//                step.entries.asSequence().map { (agentAct, occupation) ->
//                    AbstractMap.SimpleEntry(ABM.Event(time, agentAct.first, agentAct.second), occupation)
//                }
//            }
//            .flatten()
//
//    val eventVector: SparseVector<Fraction>
//    get() = events
//        .associate { (event, occupation) -> Pair(event.ordinal, Fraction(occupation)) }
//        .asVector(FractionOperators)
//
//
//
//    constructor(model: ABM<AGENT,ACT>) {
//        this.model = model
//    }
//
//    constructor(model: ABM<AGENT,ACT>, actVector: SparseVector<Fraction>): this(model) {
//        for((actId, occupation) in actVector.nonZeroEntries) {
//            this[model.eventDomain[actId]] = occupation.toDouble().roundToInt()
//        }
//    }
//
//    operator fun get(time: Int): Multiset<Pair<AGENT,ACT>> {
//        return if(time < timesteps.size) timesteps[time] else Multiset(HashMap())
//    }
//
//    inline operator fun get(time: Int, agent: AGENT, act: ACT): Int {
//        return if (time < timesteps.size) timesteps[time][Pair(agent, act)] else 0
//    }
//
//    inline operator fun get(event: ABM.Event<AGENT, ACT>): Int {
//        return get(event.time, event.agent, event.act)
//    }
//
//    inline operator fun set(time: Int, agent: AGENT, act: ACT, occupation: Int) {
//        while (time >= timesteps.size) timesteps.add(Multiset())
//        timesteps[time][Pair(agent, act)] = occupation
//    }
//
//    inline operator fun set(event: ABM.Event<AGENT, ACT>, occupation: Int) {
//        set(event.time, event.agent, event.act, occupation)
//    }
//
//
//    fun stateAt(time: Int): Multiset<AGENT> {
//        val state = Multiset<AGENT>()
//        if(time > timesteps.size) return state
//        if(time == timesteps.size) {
//            return if(timesteps.isEmpty()) state else finalState()
//        }
//        for((agentAct, occupation) in timesteps[time].entries) {
//            val (agent, _) = agentAct
//            state[agent] += occupation
//        }
//        return state
//    }
//
//
//    fun finalState(): Multiset<AGENT> {
//        val state = Multiset<AGENT>()
//        for((agentAct, occupation) in timesteps.last().entries) {
//            val (agent, act) = agentAct
//            state += model.consequences(agent, act) * occupation
//        }
//        return state
//
//    }
//
//    fun nAgents(time: Int, agent: AGENT): Int {
//        var count = 0
//        for(act in model.actDomain) {
//            count += this[time, agent, act]
//        }
//        return count
//    }
//
//
//
//    fun logPrior(): Double {
//        var logP = 0.0
//        for(t in timesteps.indices) {
//            val state = stateAt(t)
//            for((agentAct, occupation) in this[t].entries) {
//                val (agent, act) = agentAct
//                logP += occupation*ln(agent.timestep(state)[act.ordinal])
//            }
//        }
//        return logP
//    }
//
//
//
//
//    override fun toString(): String {
//        val out = StringBuilder()
//        for(step in timesteps) {
//            out.append(step.toString())
//            out.append('\n')
//        }
//        return out.toString()
//    }
//

}
