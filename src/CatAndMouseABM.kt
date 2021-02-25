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

        override fun equals(other: Any?): Boolean {
            return if(other is CatMouseAgent) {
                return other.type == type && other.position == position
            } else {
                false
            }
        }

        override fun hashCode(): Int {
            return type.hashCode()*2 + position.hashCode()
        }
    }

    class CMObservation(val agent: CatMouseAgent, val time: Int, val agentPresent: Boolean): Observation<CatMouseAgent,Acts> {

        override fun logLikelihood(trajectory: Trajectory<CatMouseAgent, Acts>): Double {
            return if(
                (agentPresent && trajectory.nAgents(time, agent) >= 1) ||
                (!agentPresent && trajectory.nAgents(time, agent) == 0)
            ) 0.0 else Double.NEGATIVE_INFINITY
        }

        override fun eventConstraints(): List<Constraint<Fraction>> {
            return if(agentPresent) {
                listOf(Constraint(hashMapOf(agent.ordinal to Fraction.ONE), ">=", Fraction.ONE).stateToEventConstraint(time))
            } else {
                listOf(Constraint(hashMapOf(agent.ordinal to Fraction.ONE), "==", Fraction.ZERO).stateToEventConstraint(time))
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

    override fun timestepEventConstraints(event: ABM.Event<CatMouseAgent, Acts>): List<Constraint<Fraction>> {
        return when(event.agent.type) {
            AgentType.CAT   -> emptyList()
            AgentType.MOUSE -> when(event.act) {
                Acts.MOVE -> listOf(
                    Constraint(
                        hashMapOf(CatMouseAgent(AgentType.CAT, event.agent.position).ordinal to Fraction.ONE),
                        ">=",
                        Fraction.ONE
                    ).stateToEventConstraint(event.time)
                )
                Acts.STAYPUT -> listOf(
                    Constraint(
                        hashMapOf(CatMouseAgent(AgentType.CAT, event.agent.position).ordinal to Fraction.ONE),
                        "==",
                        Fraction.ZERO
                    ).stateToEventConstraint(event.time)

                )
            }
        }
    }

}