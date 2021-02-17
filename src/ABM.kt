import lib.vector.SparseVector
import org.apache.commons.math3.fraction.Fraction

interface ABM<AGENT: Agent<AGENT,ACT>,ACT> {
    val actDomain: CountableDomain<ACT>
    val agentDomain: CountableDomain<AGENT>

    fun isFermionic(): Boolean = true
    fun trajectoryLogProb(x: SparseVector<Fraction>): Double = 0.0
    fun action(startState: AGENT, act: ACT): Map<AGENT,Int>
    fun timestepSupport(agent: AGENT, act: ACT): List<Constraint<Fraction>>
//    fun AGENT.timestep(otherAgents: Map<AGENT,Int>): ACT
}