object TinyABM: ABM<TinyABM.Agent,TinyABM.Acts> {

    class Agent(val x: Int)

    enum class Acts {
        MOVE,
        STAYPUT
    }

    override val actDomain = countableDomainOf<Acts>()
    override val agentDomain = object: CountableDomain<Agent> {
        override val size = 2
        override fun toIndex(agent: Agent) = agent.x
        override fun toObject(index: Int) = Agent(index)
    }

    override fun action(startState: Agent, act: Acts): Map<Agent, Int> {
        return when(act) {
            Acts.STAYPUT    -> mapOf(startState to 1)
            Acts.MOVE       -> mapOf(Agent(startState.x xor 1) to 1)
        }
    }

    override fun Agent.timestep(otherAgents: Map<Agent, Int>): Acts {
        TODO("Not yet implemented")
    }


}