import lib.collections.Multiset
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

    enum class Acts: Ordered<Acts> {
        MOVE,
        STAYPUT;

        override val domain: CountableDomain<Acts>
            get() = enumCountableDomain()
    }

    class CatMouseAgent(val type: AgentType, val position: AgentPosition): Agent<CatMouseAgent> {
        override fun timestep(others: Multiset<CatMouseAgent>): Array<Double> {
            val pCatMove = 0.25

            return if(type == AgentType.CAT) {
                arrayOf(pCatMove, 1.0-pCatMove)
            } else {
                if(others[CatMouseAgent(AgentType.CAT, position)] >= 1) {
                    arrayOf(1.0,0.0)
                } else {
                    arrayOf(0.0,1.0)
                }
            }
        }

        override val ordinal:Int
            get() = type.ordinal*2 + position.ordinal

        override val domain: CountableDomain<CatMouseAgent>
            get() = agentDomain

        override fun toString(): String {
            return "${type}:${position}"
        }
    }

    class CMObservation(val agent: CatMouseAgent, val time: Int, val agentPresent: Boolean): Observation<CatMouseAgent,Acts> {

        override fun logLikelihood(trajectory: Trajectory<CatMouseAgent, Acts>): Double {
            return if(
                (agentPresent && trajectory.nAgents(time, agent) >= 1) ||
                (!agentPresent && trajectory.nAgents(time, agent) == 0)
            ) 0.0 else Double.NEGATIVE_INFINITY
        }

        override fun constraints(): List<Constraint<Fraction>> {
            return if(agentPresent) {
                listOf(Constraint(hashMapOf(agent.ordinal to Fraction.ONE), ">=", Fraction.ONE))
            } else {
                listOf(Constraint(hashMapOf(agent.ordinal to Fraction.ONE), "==", Fraction.ZERO))
            }
        }
    }




    override val actDomain = enumCountableDomain<Acts>()
    override val agentDomain = object: CountableDomain<CatMouseAgent> {
        override val size = 4

        override fun get(index: Int) = CatMouseAgent(
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
                        hashMapOf(CatMouseAgent(AgentType.CAT, agent.position).ordinal to Fraction.ONE),
                        ">=",
                        Fraction.ONE
                    )
                )
                Acts.STAYPUT -> emptyList()
            }
        }
    }

}