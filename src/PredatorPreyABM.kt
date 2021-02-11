import org.apache.commons.math3.fraction.Fraction

object PredatorPreyABM: ABM<PredatorPreyABM.Agent, PredatorPreyABM.Acts> {
    class Agent(val x: Int, val y: Int, val isPredator: Boolean) {
        fun timestep(others: Map<Agent,Int>): Acts {
            TODO()
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

    var gridSize = 32

    override val actDomain = countableDomainOf<Acts>()

    override val agentDomain = object: CountableDomain<Agent> {
        override val size = 2 * gridSize * gridSize

        override fun toIndex(agent: Agent) =
            agent.x + agent.y * gridSize + if(agent.isPredator) gridSize * gridSize else 0

        override fun toObject(index: Int): Agent {
            val isPredator = index < gridSize*gridSize
            val x = index.rem(gridSize)
            val y = index.rem(gridSize*gridSize) / gridSize
            return Agent(x,y,isPredator)
        }

    }

    override fun action(agent: Agent, action: Acts): Map<Agent,Int> {
        return when(action) {
            Acts.DIE       -> emptyMap()
            Acts.MOVELEFT  -> mapOf(Agent(agent.x.periodicDec(gridSize),agent.y, agent.isPredator) to 1)
            Acts.MOVERIGHT -> mapOf(Agent(agent.x.periodicInc(gridSize),agent.y, agent.isPredator) to 1)
            Acts.MOVEUP    -> mapOf(Agent(agent.x,agent.y.periodicInc(gridSize), agent.isPredator) to 1)
            Acts.MOVEDOWN  -> mapOf(Agent(agent.x,agent.y.periodicDec(gridSize), agent.isPredator) to 1)
            Acts.GIVEBIRTH -> mapOf(agent to 1, Agent(agent.x.periodicInc(gridSize),agent.y, agent.isPredator) to 1)
        }
    }

    override fun Agent.timestep(otherAgents: Map<Agent, Int>) = this.timestep(otherAgents)


    // returns a constraint in terms of state occupation numbers
    fun timestepStateConstraints(agent: Agent, act: Acts): List<Constraint<Fraction>> {
        if(agent.isPredator && act == Acts.GIVEBIRTH) {
            return listOf(Constraint(hashMapOf(
                agentDomain.toIndex(Agent(agent.x.periodicDec(gridSize),agent.y, agent.isPredator)) to Fraction.ONE,
                agentDomain.toIndex(Agent(agent.x.periodicInc(gridSize),agent.y, agent.isPredator)) to Fraction.ONE,
                agentDomain.toIndex(Agent(agent.x,agent.y.periodicDec(gridSize), agent.isPredator)) to Fraction.ONE,
                agentDomain.toIndex(Agent(agent.x,agent.y.periodicInc(gridSize), agent.isPredator)) to Fraction.ONE
            ),">=", Fraction.ONE))
        }
        return emptyList()
    }


    fun Int.periodicInc(periodicity: Int): Int = (this + 1).rem(periodicity)
    fun Int.periodicDec(periodicity: Int): Int = (this + periodicity - 1).rem(periodicity)

}