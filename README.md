# AgentBasedMCMC

This is the code to support the paper [/doc/ABMCMC.pdf](./doc/ABMCMC.pdf) which shows how to perform Markov chain Monte Carlo sampling from the posterior distribution of an Agent Based Model given a definition of the agent, a prior distribution over the boundary conditions and a set of observations.

The code assumes that the model trajectory is act-Fermionic. See the paper for more details.

The code makes use of the boost and pthread libraries, so make sure these are installed.
