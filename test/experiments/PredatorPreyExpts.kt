package experiments

import PredatorPreyAgent

class PredatorPreyExpts {

    fun validTrajectory() {
        val timesteps = 2
        val continuity = ABMConstraints.continuityConstraints(timesteps, PredatorPreyAgent, PredatorPreyAgent.Companion::predatorPreyActions)

    }
}