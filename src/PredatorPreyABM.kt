import lib.collections.Multiset
import org.apache.commons.math3.fraction.Fraction
import kotlin.math.ln
import kotlin.random.Random

object PredatorPreyABM: ABM<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts> {

    enum class AgentType {
        PREDATOR,
        PREY
    }

    enum class Acts {
        DIE,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        GIVEBIRTH
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


        fun toIndex() = agentDomain.toIndex(this)
    }

    class Observation(val lookedFor: PredPreyAgent, val time: Int) {

        fun pObserveFootprints(trajectory: ABM.Trajectory<PredPreyAgent,Acts>): Double {
            val stateAtT = trajectory[time].toABMState()
            return if(stateAtT[lookedFor] >= 1) pObserveIfPresent else 0.0
        }


        companion object {
            val pObserveIfPresent = 0.9
        }
    }


    class CompletedObservation(lookedFor: PredPreyAgent, time: Int, trajectory: ABM.Trajectory<PredPreyAgent, Acts>) {
        val observation: Observation = Observation(lookedFor, time)
        val footprintsObserved: Boolean = (Random.nextDouble() < observation.pObserveFootprints(trajectory))

        fun logProb(trajectory: ABM.Trajectory<PredPreyAgent, Acts>): Double {
            return if(footprintsObserved) {
                ln(observation.pObserveFootprints(trajectory))
            } else {
                ln(1.0-observation.pObserveFootprints(trajectory))
            }
        }

        fun constraints(): List<Constraint<Fraction>> {
            return if(footprintsObserved) {
                listOf(Constraint(hashMapOf(observation.lookedFor.toIndex() to Fraction.ONE), ">=", Fraction.ONE))
            } else {
                emptyList()
            }
        }

    }

    var gridSize = 8

    override val actDomain = countableDomainOf<Acts>()

    override val agentDomain = object: CountableDomain<PredPreyAgent> {
        override val size: Int
            get() = 2 * gridSize * gridSize

        override fun toIndex(agent: PredPreyAgent) =
            agent.x + agent.y * gridSize + if(agent.type == AgentType.PREDATOR) gridSize * gridSize else 0

        override fun toObject(index: Int): PredPreyAgent {
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
                agentDomain.toIndex(PredPreyAgent(agent.x.periodicDec(),agent.y, AgentType.PREY)) to Fraction.ONE,
                agentDomain.toIndex(PredPreyAgent(agent.x.periodicInc(),agent.y, AgentType.PREY)) to Fraction.ONE,
                agentDomain.toIndex(PredPreyAgent(agent.x,agent.y.periodicDec(), AgentType.PREY)) to Fraction.ONE,
                agentDomain.toIndex(PredPreyAgent(agent.x,agent.y.periodicInc(), AgentType.PREY)) to Fraction.ONE
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