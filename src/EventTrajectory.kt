import lib.*
import java.io.Serializable
import kotlin.math.ln
import kotlin.random.Random

class EventTrajectory<AGENT>: Serializable, ArrayList<ModelEvent<AGENT>> {

    constructor(): super()

    constructor(initialCapacity: Int): super(initialCapacity)

    fun toStateTrajectory(): List<Multiset<AGENT>> {
        return this.map { it.consequences() }
    }

    fun impliedStartState(): Multiset<AGENT> {
        return this.first().primaryRequirements()
    }


//    fun initialPrimaryRequirements(): Multiset<AGENT> {
//        if(size == 0) return emptyMultiset()
//        return first().primaryRequirements()
//    }


    fun generateObservations(pObserve: Double): List<Multiset<AGENT>> {
        return this.map { modelEvent ->
            val observedState = HashMultiset<AGENT>()
            modelEvent.consequences.filterTo(observedState) { Random.nextDouble() < pObserve }
            observedState
        }
    }


    fun logProb() = sumByDouble { it.logProb() }


    fun isFeasible(startState: Multiset<AGENT>, isPartial: Boolean): Boolean {
        var state = startState
        this.forEach { modelEvent ->
            if(!modelEvent.isSatisfiedBy(state, isPartial)) return false
            state = modelEvent.consequences
        }
        return true
    }

    companion object {
        private const val serialVersionUID: Long = -7218603259279138400
    }

}