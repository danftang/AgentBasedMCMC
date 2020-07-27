package lib

import Agent
import Predator
import Prey
import ModelEvent
import java.lang.IllegalStateException

fun Multiset<Agent>.divergence(other: Multiset<Agent>, gridSize: Int): Double {
    val predators = KDTree.Manhattan<Boolean>(2)
    val prey = KDTree.Manhattan<Boolean>(2)
    other.supportSet.forEach { agent ->
        val x = agent.pos.rem(gridSize).toDouble()
        val y = agent.pos.div(gridSize).toDouble()
        when(agent) {
            is Predator -> predators[x, y] = true
            is Prey -> prey[x, y] = true
        }
    }
    return if(size == 0) 0.0 else this.sumByDouble { agent ->
        val x = agent.pos.rem(gridSize).toDouble()
        val y = agent.pos.div(gridSize).toDouble()
        when(agent) {
            is Predator -> predators.nearestNeighbour(x, y).distance
            is Prey -> prey.nearestNeighbour(x, y).distance
            else -> throw(IllegalStateException("unrecognized agent type"))
        }
    }/this.size
}


fun Multiset<Agent>.satisfies(modelEvent: ModelEvent<Agent>, isPartial: Boolean) = modelEvent.isSatisfiedBy(this, isPartial)


fun Multiset<Agent>.comparisonPlot(other: Multiset<Agent>, gridSize: Int) {
    val thisData = this.supportSet.asSequence().flatMap { agent ->
        sequenceOf(agent.pos.rem(gridSize), agent.pos.div(gridSize), if(agent is Prey) 1 else 2)
    }
    val otherData = other.supportSet.asSequence().flatMap { agent ->
        sequenceOf(agent.pos.rem(gridSize), agent.pos.div(gridSize), if(agent is Prey) 1 else 2)
    }

    gnuplot {
        invoke("set linetype 1 lc 'red'") // prey are red
        invoke("set linetype 2 lc 'blue'") // predator are blue

        // plot this with squares and other with triangles
        val thisDoc = heredoc(thisData,3)
        val otherDoc = heredoc(otherData,3)
        invoke("plot [0:$gridSize][0:$gridSize] $thisDoc using ($1+0.5):($2+0.5):3 with points pointsize 1.8  pointtype 19 lc variable notitle")
        invoke("replot $otherDoc using ($1+0.5):($2+0.5):3 with points pointsize 1.8  pointtype 23 lc variable notitle")
    }
}



