interface ABM<AGENT,ACT> {
    val actDomain: CountableDomain<ACT>
    val agentDomain: CountableDomain<AGENT>

    fun action(startState: AGENT, act: ACT): Map<AGENT,Int>
    fun AGENT.timestep(otherAgents: Map<AGENT,Int>): ACT
}