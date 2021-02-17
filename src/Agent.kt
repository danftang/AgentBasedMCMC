interface Agent<AGENT,ACT> {
    fun timestep(others: Map<Int,AGENT>): ACT
}