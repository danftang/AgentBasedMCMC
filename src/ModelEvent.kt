import lib.HashMultiset
import lib.Multiset
import lib.isDisjoint
import java.io.Serializable
import kotlin.math.ln

class ModelEvent<AGENT>: Serializable, HashMultiset<Event<AGENT>>() {
    val consequences: Multiset<AGENT>
        get() {
            val allConsequences = HashMultiset<AGENT>()
            this.forEach {
                allConsequences.addAll(it.consequences)
            }
            return allConsequences
        }

    fun isSatisfiedBy(state: Multiset<AGENT>, isPartial: Boolean): Boolean {
        val primaryRequirements = HashMultiset<AGENT>()
        this.forEach { event ->
            if(!event.absenceRequirements.isDisjoint(state)) return false
            if(!state.containsAll(event.presenceRequirements)) return false
            primaryRequirements.add(event.primaryAgent)
        }
        if(isPartial) {
            if(!primaryRequirements.isSubsetOf(state)) return false
        } else {
            if (state != primaryRequirements) return false
        }
        return true
    }

    fun logProb() = this.sumByDouble {  event -> ln(event.rate) }
}