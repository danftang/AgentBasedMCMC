import org.apache.commons.math3.fraction.Fraction

object TinyABM: ABM<TinyABM.TinyAgent,TinyABM.Acts> {

    class TinyAgent(val x: Int): Agent<TinyAgent, Acts> {
        override fun timestep(others: Map<Int, TinyAgent>): Acts {
            TODO("Not yet implemented")
        }
    }

    enum class Acts {
        MOVE,
        STAYPUT
    }

    override val actDomain = countableDomainOf<Acts>()
    override val agentDomain = object: CountableDomain<TinyAgent> {
        override val size = 2
        override fun toIndex(agent: TinyAgent) = agent.x
        override fun toObject(index: Int) = TinyAgent(index)
    }

    override fun action(startState: TinyAgent, act: Acts): Map<TinyAgent, Int> {
        return when(act) {
            Acts.STAYPUT    -> mapOf(startState to 1)
            Acts.MOVE       -> mapOf(TinyAgent(startState.x xor 1) to 1)
        }
    }

    override fun timestepSupport(agent: TinyAgent, act: Acts): List<Constraint<Fraction>> {
        return emptyList()
    }
}