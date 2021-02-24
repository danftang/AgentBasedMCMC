import lib.collections.Multiset
import org.apache.commons.math3.fraction.Fraction
import kotlin.math.ln
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
            val pPredBirthGivenPrey = 0.1 // birth prob given prey
            val pPredDie = 0.1 // death prob
            val pPreyBirth = 0.1 // birth prob
            val pPreyDie = 0.1 // death prob
            val pPreyEatenGivenPred = 0.2 // death prob given pred

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

        override fun constraints(): List<Constraint<Fraction>> {
            return if(footprintsObserved) {
                listOf(Constraint(hashMapOf(lookedFor.ordinal to Fraction.ONE), ">=", Fraction.ONE))
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

    override fun action(startState: PredPreyAgent, act: Acts): Map<PredPreyAgent,Int> {
        return when(act) {
            Acts.DIE       -> emptyMap()
            Acts.MOVELEFT  -> mapOf(PredPreyAgent(startState.x.periodicDec(),startState.y, startState.type) to 1)
            Acts.MOVERIGHT -> mapOf(PredPreyAgent(startState.x.periodicInc(),startState.y, startState.type) to 1)
            Acts.MOVEUP    -> mapOf(PredPreyAgent(startState.x,startState.y.periodicInc(), startState.type) to 1)
            Acts.MOVEDOWN  -> mapOf(PredPreyAgent(startState.x,startState.y.periodicDec(), startState.type) to 1)
            Acts.GIVEBIRTH -> mapOf(startState to 1, PredPreyAgent(startState.x.periodicInc(),startState.y, startState.type) to 1)
        }
    }


    // returns a constraint in terms of state occupation numbers
    // Manually generated for now...
    override fun timestepSupport(agent: PredPreyAgent, act: Acts): List<Constraint<Fraction>> {
        if(agent.type == AgentType.PREDATOR && act == Acts.GIVEBIRTH) {
            return listOf(Constraint(hashMapOf(
                PredPreyAgent(agent.x.periodicDec(),agent.y, AgentType.PREY).ordinal to Fraction.ONE,
                PredPreyAgent(agent.x.periodicInc(),agent.y, AgentType.PREY).ordinal to Fraction.ONE,
                PredPreyAgent(agent.x,agent.y.periodicDec(), AgentType.PREY).ordinal to Fraction.ONE,
                PredPreyAgent(agent.x,agent.y.periodicInc(), AgentType.PREY).ordinal to Fraction.ONE
            ),">=", Fraction.ONE))
        }
        return emptyList()
    }


    fun Int.periodicInc(): Int = (this + 1).rem(gridSize)
    fun Int.periodicDec(): Int = (this + gridSize - 1).rem(gridSize)

    fun randomState(pPredator: Double, pPrey: Double): Multiset<PredPreyAgent> {
        val state = Multiset<PredPreyAgent>()
        for(x in 0 until gridSize) {
            for(y in 0 until gridSize) {
                if(Random.nextDouble() < pPredator) state[PredPreyAgent(x,y,AgentType.PREDATOR)] = 1
                if(Random.nextDouble() < pPrey)     state[PredPreyAgent(x,y,AgentType.PREY)] = 1
            }
        }
        return state
    }


}