interface ABMDomain<AGENT,ACT: Enum<ACT>> {
    val nAgentStates: Int

    fun toIndex(agent: AGENT): Int
    fun toAgent(index: Int): AGENT
    fun actValues(): Array<ACT>
}