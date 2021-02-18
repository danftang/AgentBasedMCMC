import org.apache.commons.math3.fraction.Fraction

object CatAndMouseABM: ABM<CatAndMouseABM.CatMouseAgent,CatAndMouseABM.Acts> {

    enum class AgentType {
        CAT,
        MOUSE
    }

    enum class AgentPosition {
        LEFT,
        RIGHT;

        fun notHere() = if(this == LEFT) RIGHT else LEFT
    }

    enum class Acts {
        MOVE,
        STAYPUT
    }

    class CatMouseAgent(val type: AgentType, val position: AgentPosition): Agent<CatMouseAgent, Acts> {

        override fun timestep(others: Map<Int, CatMouseAgent>): Acts {
            TODO("Not yet implemented")
        }

        fun toIndex() = agentDomain.toIndex(this)
    }



    override val actDomain = countableDomainOf<Acts>()
    override val agentDomain = object: CountableDomain<CatMouseAgent> {
        override val size = 4
        override fun toIndex(agent: CatMouseAgent) = agent.type.ordinal*2 + agent.position.ordinal
        override fun toObject(index: Int) = CatMouseAgent(
            AgentType.values()[index.shr(1)],
            AgentPosition.values()[index and 1]
        )
    }

    override fun action(startState: CatMouseAgent, act: Acts): Map<CatMouseAgent, Int> {
        return when(act) {
            Acts.STAYPUT    -> mapOf(startState to 1)
            Acts.MOVE       -> mapOf(CatMouseAgent(startState.type, startState.position.notHere()) to 1)
        }
    }

    override fun timestepSupport(agent: CatMouseAgent, act: Acts): List<Constraint<Fraction>> {
        return when(agent.type) {
            AgentType.CAT   -> emptyList()
            AgentType.MOUSE -> when(act) {
                Acts.MOVE -> listOf(
                    Constraint(
                        hashMapOf(CatMouseAgent(AgentType.CAT, agent.position).toIndex() to Fraction.ONE),
                        ">=",
                        Fraction.ONE
                    )
                )
                Acts.STAYPUT -> emptyList()
            }
        }
    }
}