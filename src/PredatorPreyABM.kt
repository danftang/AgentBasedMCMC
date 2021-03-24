import lib.Gnuplot
import lib.collections.Multiset
import lib.collections.multisetOf
import lib.gnuplot
import org.apache.commons.math3.fraction.Fraction
import kotlin.math.ln
import kotlin.math.roundToInt
import kotlin.random.Random

object PredatorPreyABM: ABM<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts> {

    enum class AgentType {
        PREDATOR,
        PREY
    }

    enum class Acts: Ordered<Acts> {
        DIE,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        GIVEBIRTH;

        override val domain: CountableDomain<Acts>
            get() = enumCountableDomain()
    }


    class PredPreyAgent(val x: Int, val y: Int, val type: AgentType): Agent<PredPreyAgent> {

        override fun timestep(others: Multiset<PredPreyAgent>): Array<Double> {
            val pPredBirthGivenPrey = 0.5 // birth prob given prey
            val pPredDie = 0.07 // death prob
            val pPreyBirth = 0.06 // birth prob
            val pPreyDie = 0.03 // death prob
            val pPreyEatenGivenPred = 0.55 // death prob given pred

            val actDistribution = Array(Acts.values().size) { 0.0 }
            if(type == AgentType.PREDATOR) {
                actDistribution[Acts.DIE.ordinal] = pPredDie
                if(surroundingCountOf(AgentType.PREY, others) >= 1) {
                    actDistribution[Acts.GIVEBIRTH.ordinal] = pPredBirthGivenPrey
                }
            } else {
                actDistribution[Acts.GIVEBIRTH.ordinal] = pPreyBirth
                actDistribution[Acts.DIE.ordinal] = pPreyDie
                if(surroundingCountOf(AgentType.PREDATOR, others) >= 1) {
                    actDistribution[Acts.DIE.ordinal] += pPreyEatenGivenPred
                }
            }
            val moveProb = 0.25*(1.0 - actDistribution[Acts.DIE.ordinal] - actDistribution[Acts.GIVEBIRTH.ordinal])
            actDistribution[Acts.MOVEUP.ordinal] = moveProb
            actDistribution[Acts.MOVEDOWN.ordinal] = moveProb
            actDistribution[Acts.MOVELEFT.ordinal] = moveProb
            actDistribution[Acts.MOVERIGHT.ordinal] = moveProb
            return actDistribution
        }


        fun surroundingCountOf(type: AgentType, others: Multiset<PredPreyAgent>): Int =
            others[PredPreyAgent(x.periodicInc(),y,type)] +
                    others[PredPreyAgent(x.periodicDec(),y,type)] +
                    others[PredPreyAgent(x,y.periodicInc(),type)] +
                    others[PredPreyAgent(x,y.periodicDec(),type)]


        override val ordinal: Int
            get() = x + y * gridSize + if(type == AgentType.PREDATOR) gridSize * gridSize else 0

        override val domain: CountableDomain<PredPreyAgent>
            get() = agentDomain

        override fun equals(other: Any?): Boolean {
            return if(other is PredPreyAgent) {
                x == other.x && y == other.y && type == other.type
            } else false
        }

        override fun hashCode() = ordinal

        override fun toString(): String {
            return "$type($x,$y)"
        }

    }


    class PPObservation(val time: Int, val lookedFor: PredPreyAgent, trajectory: Trajectory<PredPreyAgent, Acts>): Observation<PredPreyAgent, Acts> {
        val footprintsObserved: Boolean = (Random.nextDouble() < pObserveFootprints(time, lookedFor, trajectory))

        override fun logLikelihood(trajectory: Trajectory<PredPreyAgent, Acts>): Double {
            return if(footprintsObserved) {
                ln(pObserveFootprints(time, lookedFor, trajectory))
            } else {
                ln(1.0-pObserveFootprints(time, lookedFor, trajectory))
            }
        }

        override fun eventConstraints(): List<MutableConstraint<Fraction>> {
            return if(footprintsObserved) {
                listOf(MutableConstraint(hashMapOf(lookedFor.ordinal to Fraction.ONE), ">=", Fraction.ONE).stateToEventConstraint(time))
            } else {
                emptyList()
            }
        }

        companion object {
            val pObserveIfPresent = 0.9

            fun pObserveFootprints(time: Int, lookedFor: PredPreyAgent, trajectory: Trajectory<PredPreyAgent, Acts>): Double {
                return if(trajectory.nAgents(time, lookedFor) >= 1) pObserveIfPresent else 0.0
            }
        }

    }


    // Prior knowledge that the occupation of each state at time t=0 was drawn from a Binomial, expressed as an observation.
    class Prior(val initialPredatorDensity: Double, val initialPreyDensity: Double): Observation<PredPreyAgent, Acts> {
        override fun logLikelihood(trajectory: Trajectory<PredPreyAgent, Acts>): Double {
            var logP = 0.0
            var binomialLogP: Pair<Double,Double>
            val predLogP = Pair(ln(initialPredatorDensity),ln(1.0-initialPredatorDensity))
            val preyLogP = Pair(ln(initialPreyDensity),ln(1.0-initialPreyDensity))
            for(agent in agentDomain) {
                binomialLogP = if(agent.type == AgentType.PREDATOR) predLogP else preyLogP
                logP += if(trajectory.nAgents(0,agent) > 0) binomialLogP.first else binomialLogP.second
            }
            return logP
        }

        override fun eventConstraints(): List<MutableConstraint<Fraction>> = emptyList() // no states are impossible

    }


    var gridSize = 8

    override val actDomain = enumCountableDomain<Acts>()

    override val agentDomain = object: CountableDomain<PredPreyAgent> {
        override val size: Int
            get() = 2 * gridSize * gridSize

        override fun get(index: Int): PredPreyAgent {
            val type = if(index < gridSize*gridSize) AgentType.PREY else AgentType.PREDATOR
            val x = index.rem(gridSize)
            val y = index.rem(gridSize*gridSize) / gridSize
            return PredPreyAgent(x,y,type)
        }

    }

    override fun consequences(startState: PredPreyAgent, act: Acts): Multiset<PredPreyAgent> {
        return when(act) {
            Acts.DIE       -> multisetOf()
            Acts.MOVELEFT  -> multisetOf(PredPreyAgent(startState.x.periodicDec(),startState.y, startState.type) to 1)
            Acts.MOVERIGHT -> multisetOf(PredPreyAgent(startState.x.periodicInc(),startState.y, startState.type) to 1)
            Acts.MOVEUP    -> multisetOf(PredPreyAgent(startState.x,startState.y.periodicInc(), startState.type) to 1)
            Acts.MOVEDOWN  -> multisetOf(PredPreyAgent(startState.x,startState.y.periodicDec(), startState.type) to 1)
            Acts.GIVEBIRTH -> multisetOf(startState to 1, PredPreyAgent(startState.x.periodicInc(),startState.y, startState.type) to 1)
        }
    }


    // returns a constraint in terms of state occupation numbers
    // Manually generated for now...
    override fun timestepEventConstraints(event: ABM.Event<PredPreyAgent, Acts>): List<MutableConstraint<Fraction>> {
        if(event.agent.type == AgentType.PREDATOR && event.act == Acts.GIVEBIRTH) {
            return listOf(MutableConstraint(hashMapOf(
                PredPreyAgent(event.agent.x.periodicDec(),event.agent.y, AgentType.PREY).ordinal to Fraction.ONE,
                PredPreyAgent(event.agent.x.periodicInc(),event.agent.y, AgentType.PREY).ordinal to Fraction.ONE,
                PredPreyAgent(event.agent.x,event.agent.y.periodicDec(), AgentType.PREY).ordinal to Fraction.ONE,
                PredPreyAgent(event.agent.x,event.agent.y.periodicInc(), AgentType.PREY).ordinal to Fraction.ONE
            ),">=", Fraction.ONE).stateToEventConstraint(event.time))
        }
        return emptyList()
    }


    fun Int.periodicInc(): Int = (this + 1).rem(gridSize)
    fun Int.periodicDec(): Int = (this + gridSize - 1).rem(gridSize)


    // Occupation of each state is drawn from a Binomial distribution
    fun randomFermionicState(pPredator: Double, pPrey: Double): Multiset<PredPreyAgent> {
        val state = Multiset<PredPreyAgent>()
        for(x in 0 until gridSize) {
            for(y in 0 until gridSize) {
                if(Random.nextDouble() < pPredator) state[PredPreyAgent(x,y,AgentType.PREDATOR)] = 1
                if(Random.nextDouble() < pPrey)     state[PredPreyAgent(x,y,AgentType.PREY)] = 1
            }
        }
        return state
    }


    fun plotHeatMap(state: Multiset<PredPreyAgent>): Gnuplot {
        val maxCount = state.occupationNumbers.max()?:1

        val rabbitData = agentDomain.asSequence()
            .filter { it.type == AgentType.PREY }
            .flatMap { agent ->
//                sequenceOf<Number>(agent.x, agent.y, (state[agent]*255.0/maxCount).roundToInt(), 0, 0, 128)
                sequenceOf<Number>(agent.x, agent.y, (ln(state[agent]+1.0)*255.0/ln(maxCount+1.0)).roundToInt(), 0, 0, 128)
            }

        val foxData = agentDomain.asSequence()
            .filter { it.type == AgentType.PREDATOR }
            .flatMap { agent ->
//                sequenceOf<Number>(agent.x, agent.y, 0, 0, (state[agent]*255.0/maxCount).roundToInt(), 128)
                sequenceOf<Number>(agent.x, agent.y, 0, 0, (ln(state[agent]+1.0)*255.0/ln(maxCount+1.0)).roundToInt(), 128)
            }

//        val rabbitData =  state.entries.asSequence()
//            .filter { it.key.type == PredatorPreyABM.AgentType.PREY }
//            .flatMap { (agent, count) ->
//                sequenceOf<Number>(agent.x, agent.y, (count*255.0/maxCount).roundToInt(), 0, 0, 128)
//            }
//
//        val foxData = state.entries.asSequence()
//            .filter { it.key.type == PredatorPreyABM.AgentType.PREDATOR }
//            .flatMap { (agent, count) ->
//                sequenceOf<Number>(agent.x, agent.y, 0, 0, (count*255.0/maxCount).roundToInt(), 128)
//            }

        val gp = Gnuplot()
        gnuplot {  }
        with(gp) {
            val rData = heredoc(rabbitData,6)
            val fData = heredoc(foxData,6)
            invoke("plot [-0.5:${gridSize-0.5}][-0.5:${gridSize-0.5}] $rData with rgbalpha")
            invoke("replot $fData with rgbalpha")
        }
        return gp
    }


    fun Gnuplot.replotPoints(state: Multiset<PredPreyAgent>): Gnuplot {
        val stateData = state.entries.asSequence().flatMap { (agent, _) ->
            sequenceOf(agent.x, agent.y, if (agent.type == AgentType.PREY) 1 else 2)
        }.toList()
        invoke("set linetype 1 lc 'red'")
        invoke("set linetype 2 lc 'blue'")
        if (stateData.isNotEmpty()) {
            val pointData = heredoc(stateData, 3)
            invoke("replot $pointData with points pointtype 5 pointsize 0.5 lc variable")
        }
        return this
    }


}