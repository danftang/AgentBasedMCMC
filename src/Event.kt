import lib.Multiset
import lib.MutableMultiset
import lib.toMutableMultiset
import java.io.Serializable
import java.lang.IllegalArgumentException

class Event<AGENT>(val presenceRequirements: Set<AGENT>, val absenceRequirements: Set<AGENT>, val consequences: Multiset<AGENT>, val rate: Double, val primaryAgent: AGENT): Serializable {
    val additions: MutableMultiset<AGENT>
        get() = consequences - presenceRequirements
    val secondaryRequirements: Set<AGENT>
        get() = presenceRequirements.minus(primaryAgent)
//    val totalRequirements: Set<AGENT>
//        get() = requirements.union(absenceRequirements)

//    val deletions: Set<AGENT>
//        get() = requirements - consequences.counts.keys


    fun delta(agent: AGENT): Int {
        return consequences.count(agent) - (if(presenceRequirements.contains(agent)) 1 else 0)
    }

    fun rateFor(state: Multiset<AGENT>): Double {
        return rateFor(state.supportSet)
//        return requirements.fold(rate) {r, agent -> r*state.count(agent)}
    }

    fun rateFor(state: Set<AGENT>): Double {
        absenceRequirements.forEach { if(state.contains(it)) return 0.0 }
        if(!state.containsAll(presenceRequirements)) return 0.0
        return rate
    }


    fun actOn(state: Multiset<AGENT>): MutableMultiset<AGENT> {
        val result = state.toMutableMultiset()
        if(!result.contains(primaryAgent)) throw(IllegalArgumentException("Applying act to state without primary agent"))
        deltas().forEach {entry ->
            if(entry.value > 0) {
                result.add(entry.key, entry.value)
            } else {
                if(!result.remove(entry.key, -entry.value))
                    throw(IllegalArgumentException("Applying act to state without enough agents"))
            }
        }
        return result
    }

    fun modifiedAgents(): Set<AGENT> {
        return additions.supportSet.union(setOf(primaryAgent))
    }

    fun deltas(): Map<AGENT,Int> {
        val d = HashMap<AGENT,Int>()
        presenceRequirements.forEach {
            d[it] = -1
        }
        consequences.counts.forEach {
            d.compute(it.key) { _, oldCount ->
                val newCount = (oldCount?:0) + it.value
                if(newCount == 0) null else newCount
            }
        }
        return d
    }

    override fun toString(): String {
        return "$primaryAgent${secondaryRequirements}${absenceRequirements} -> ${consequences}"
    }
}