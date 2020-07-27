import lib.toMutableMultiset
import java.io.Serializable
import kotlin.math.abs

abstract class Agent(val pos: Int): Serializable {
    abstract fun copyAt(id: Int): Agent




    fun die(h: Hamiltonian<Agent>, rate: Double) {
        h += action(rate)
    }

    fun right(size: Int) = pos - pos.rem(size) + (pos+1).rem(size)
    fun left(size: Int) = pos - pos.rem(size) + (pos+size-1).rem(size)
    fun up(size: Int) = (pos + size).rem(size*size)
    fun down(size: Int) = (pos + size*size - size).rem(size*size)


    fun action(rate: Double, vararg addedAgents: Agent) : Event<Agent> {
        return Event(setOf(this), emptySet(), addedAgents.toMutableMultiset(), rate, this)
    }


    fun action(rate: Double, absences: Set<Agent>, vararg addedAgents: Agent) : Event<Agent> {
        return Event(setOf(this), absences, addedAgents.toMutableMultiset(), rate, this)
    }


    fun interaction(rate: Double, otherAgent: Agent, vararg addedAgents: Agent) : Event<Agent> {
        val added = addedAgents.toMutableMultiset()
        added.add(otherAgent) // require other agent is unchanged
        return Event(setOf(this, otherAgent), emptySet(), added, rate, this)
    }


    fun interaction(rate: Double, otherAgent: Agent, absences: Set<Agent>, vararg addedAgents: Agent) : Event<Agent> {
        val added = addedAgents.toMutableMultiset()
        added.add(otherAgent) // require other agent is unchanged
        return Event(setOf(this, otherAgent), absences, added, rate, this)
    }

//    fun translate(transVector: Int, gridSize: Int): Agent {
//        val newX = (pos + transVector + gridSize).rem(gridSize)
//        val newY = (transVector.div(gridSize) + pos.div(gridSize) + gridSize).rem(gridSize)
//        return copyAt(newY*gridSize + newX)
//    }
}