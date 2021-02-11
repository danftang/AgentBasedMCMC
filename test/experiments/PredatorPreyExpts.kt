package experiments

import Constraint
import PredatorPreyABM
import org.apache.commons.math3.fraction.Fraction

class PredatorPreyExpts {

    fun fermionicPredPrey() {
        val timesteps = 2
        val continuity = ABMConstraints.continuityConstraints(timesteps, PredatorPreyABM)
        val fermionic = ABMConstraints.fermionicConstraints(timesteps, PredatorPreyABM)

        val agentConstraints = ArrayList<Constraint<Fraction>>()
        for(state in 0 until PredatorPreyABM.agentDomain.size) {
            for(act in 0 until PredatorPreyABM.actDomain.size) {
                val timestepConstraints = PredatorPreyABM.timestepStateConstraints(
                    PredatorPreyABM.agentDomain.toObject(state),
                    PredatorPreyABM.actDomain.toObject(act)
                )
                for(t in 0 until timesteps) {
                    agentConstraints.addAll(timestepConstraints
                        .map { ABMConstraints.stateConstraintToActConstraint(it, t, PredatorPreyABM) }
                        .map { ABMConstraints.fermionicXImpliesY(state, it) }
                    )
                }
            }
        }

    }
}