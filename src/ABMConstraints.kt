import lib.abstractAlgebra.IntOperators
import lib.sparseMatrix.GridMapMatrix

object ABMConstraints {
    // nStates      - number of agent states
    // nTimesteps   - number of timesteps
    // nActs        - number of agent actions
    // action       - The "action function", takes an agent state, psi, and an action, a, and returns a sparse vector X
    //                such that x_i is the number of agents in state i after an agent in state psi performs act a.
    fun<ACT: Enum<ACT>,AGENT> continuityConstraints(nTimesteps: Int, domain: ABMDomain<AGENT,ACT>, action: (AGENT, ACT) -> Map<AGENT,Int>): GridMapMatrix<Int> {
        val nStates = domain.nAgentStates
        val acts = domain.actValues()
        val nActs = acts.size
        val constraints = GridMapMatrix(
            IntOperators,
            (nTimesteps-1)*nStates,
            nTimesteps*nStates*nActs
        )

        // first do leaving edges
        var i = 0
        for(t in 1 until nTimesteps) {
            for(state in 0 until nStates) {
                for(act in acts) {
                    constraints[i, t*nStates*nActs + state*nActs + act.ordinal] = -1
                }
                i++
            }
        }
        // now do incoming edges
        for(state in 0 until nStates) {
            for(act in acts) {
                val consequences = action(domain.toAgent(state), act)
                for((resultState, n) in consequences) {
                    for(t in 0 until nTimesteps-1) {
                        constraints[t*nStates + domain.toIndex(resultState), t*nStates*nActs + state*nActs + act.ordinal] = n
                    }
                }
            }
        }
        return constraints
    }


    fun fermionicConstraints(): GridMapMatrix<Int> {

    }

}