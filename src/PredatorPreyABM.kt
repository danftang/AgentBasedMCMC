import org.apache.commons.math3.fraction.Fraction

object PredatorPreyABM: ABM<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts> {
    class PredPreyAgent(val x: Int, val y: Int, val isPredator: Boolean): Agent<PredPreyAgent,Acts> {
        override fun timestep(others: Map<Int, PredPreyAgent>): Acts {
            TODO("Not yet implemented")
        }
    }

    enum class Acts {
        DIE,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        GIVEBIRTH
    }

    var gridSize = 8

    override val actDomain = countableDomainOf<Acts>()

    override val agentDomain = object: CountableDomain<PredPreyAgent> {
        override val size: Int
            get() = 2 * gridSize * gridSize

        override fun toIndex(agent: PredPreyAgent) =
            agent.x + agent.y * gridSize + if(agent.isPredator) gridSize * gridSize else 0

        override fun toObject(index: Int): PredPreyAgent {
            val isPredator = index < gridSize*gridSize
            val x = index.rem(gridSize)
            val y = index.rem(gridSize*gridSize) / gridSize
            return PredPreyAgent(x,y,isPredator)
        }

    }

    override fun action(agent: PredPreyAgent, action: Acts): Map<PredPreyAgent,Int> {
        return when(action) {
            Acts.DIE       -> emptyMap()
            Acts.MOVELEFT  -> mapOf(PredPreyAgent(agent.x.periodicDec(gridSize),agent.y, agent.isPredator) to 1)
            Acts.MOVERIGHT -> mapOf(PredPreyAgent(agent.x.periodicInc(gridSize),agent.y, agent.isPredator) to 1)
            Acts.MOVEUP    -> mapOf(PredPreyAgent(agent.x,agent.y.periodicInc(gridSize), agent.isPredator) to 1)
            Acts.MOVEDOWN  -> mapOf(PredPreyAgent(agent.x,agent.y.periodicDec(gridSize), agent.isPredator) to 1)
            Acts.GIVEBIRTH -> mapOf(agent to 1, PredPreyAgent(agent.x.periodicInc(gridSize),agent.y, agent.isPredator) to 1)
        }
    }


    // returns a constraint in terms of state occupation numbers
    override fun timestepSupport(agent: PredPreyAgent, act: Acts): List<Constraint<Fraction>> {
        if(agent.isPredator && act == Acts.GIVEBIRTH) {
            return listOf(Constraint(hashMapOf(
                agentDomain.toIndex(PredPreyAgent(agent.x.periodicDec(gridSize),agent.y, agent.isPredator)) to Fraction.ONE,
                agentDomain.toIndex(PredPreyAgent(agent.x.periodicInc(gridSize),agent.y, agent.isPredator)) to Fraction.ONE,
                agentDomain.toIndex(PredPreyAgent(agent.x,agent.y.periodicDec(gridSize), agent.isPredator)) to Fraction.ONE,
                agentDomain.toIndex(PredPreyAgent(agent.x,agent.y.periodicInc(gridSize), agent.isPredator)) to Fraction.ONE
            ),">=", Fraction.ONE))
        }
        return emptyList()
    }


    fun Int.periodicInc(periodicity: Int): Int = (this + 1).rem(periodicity)
    fun Int.periodicDec(periodicity: Int): Int = (this + periodicity - 1).rem(periodicity)

}