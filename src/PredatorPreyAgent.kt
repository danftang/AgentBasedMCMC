class PredatorPreyAgent(val x: Int, val y: Int, val isPredator: Boolean) {
    enum class ACTDOMAIN {
        DIE,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        GIVEBIRTH
    }


    companion object: ABMDomain<PredatorPreyAgent, ACTDOMAIN> {
        var gridSize = 32
        override val nAgentStates: Int
            get() = 2 * gridSize * gridSize


        override fun toIndex(agent: PredatorPreyAgent): Int {
            return agent.x + agent.y * gridSize + if(agent.isPredator) gridSize * gridSize else 0
        }


        override fun toAgent(index: Int): PredatorPreyAgent {
            val isPredator = index < gridSize*gridSize
            val x = index.rem(gridSize)
            val y = index.rem(gridSize*gridSize) / gridSize
            return PredatorPreyAgent(x,y,isPredator)
        }


        override fun actValues() = ACTDOMAIN.values()


        fun predatorPreyActions(agent: PredatorPreyAgent, action: ACTDOMAIN): Map<PredatorPreyAgent,Int> {
            when(action) {
                ACTDOMAIN.DIE       -> return emptyMap()
                ACTDOMAIN.MOVELEFT  -> return mapOf(PredatorPreyAgent((agent.x+gridSize-1).rem(gridSize),agent.y, agent.isPredator) to 1)
                ACTDOMAIN.MOVERIGHT -> return mapOf(PredatorPreyAgent((agent.x+1).rem(gridSize),agent.y, agent.isPredator) to 1)
                ACTDOMAIN.MOVEUP    -> return mapOf(PredatorPreyAgent(agent.x,(agent.y+1).rem(gridSize), agent.isPredator) to 1)
                ACTDOMAIN.MOVEDOWN  -> return mapOf(PredatorPreyAgent(agent.x,(agent.y+gridSize-1).rem(gridSize), agent.isPredator) to 1)
                ACTDOMAIN.GIVEBIRTH -> return mapOf(agent to 1, PredatorPreyAgent((agent.x+1).rem(gridSize),agent.y, agent.isPredator) to 1)
            }
        }

    }
}