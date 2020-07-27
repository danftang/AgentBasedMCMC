import java.io.Serializable

open class Hamiltonian<AGENT>: ArrayList<Event<AGENT>>(), Serializable {
    val requirementIndex: Map<AGENT, Set<Event<AGENT>>>          by lazy { indexBy { it.presenceRequirements } }
    val primaryRequirementIndex: Map<AGENT, Set<Event<AGENT>>>   by lazy { indexBy { setOf(it.primaryAgent) } }
    val consequenceIndex: Map<AGENT, Set<Event<AGENT>>>          by lazy { indexBy { it.consequences } }
    val secondaryRequirementIdnex: Map<AGENT, Set<Event<AGENT>>> by lazy { indexBy { it.presenceRequirements.minus(it.primaryAgent) } }
    val deltaIndex: Map<AGENT, Set<Event<AGENT>>>                by lazy { indexBy {
            it.presenceRequirements.minus(it.consequences).union(it.consequences.minus(it.presenceRequirements))
        }
    }
    val allStates: Set<AGENT>
        get() = primaryRequirementIndex.keys

    fun eventsWithRequirement(agent: AGENT) = requirementIndex[agent]?:emptySet()
    fun eventsWithPrimaryRequirement(agent: AGENT) = primaryRequirementIndex[agent]?:emptySet()
    fun eventsWithSecondaryRequirement(agent: AGENT) = secondaryRequirementIdnex[agent]?:emptySet()
    fun eventsWithConsequence(agent: AGENT) = consequenceIndex[agent]?:emptySet()
    fun eventsThatChange(agent: AGENT) = deltaIndex[agent]?:emptySet()


    fun eventsPresenceSatisfiedBy(footprint: Set<AGENT>): List<Event<AGENT>> {
        return footprint.flatMap { agent ->
            eventsWithPrimaryRequirement(agent)
                .filter { event -> footprint.containsAll(event.presenceRequirements) }

        }
    }


    fun <K>indexBy(keySetSelector: (Event<AGENT>) -> Set<K>): Map<K,Set<Event<AGENT>>> {
        val index = HashMap<K, MutableSet<Event<AGENT>>>()
        this.forEach { event ->
            keySetSelector(event).forEach { key ->
                index.getOrPut(key,{HashSet()}).add(event)
            }
        }
        return index
    }
}