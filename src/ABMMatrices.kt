import lib.abstractAlgebra.IntOperators
import lib.sparseIntMatrix.HashColIntMatrix
import lib.sparseIntMatrix.SparseIntVector
import lib.gnuplot
import lib.sparseMatrix.GridMapMatrix

object ABMMatrices {
    fun twoDabmMatrix(gridSize: Int, timesteps: Int): GridMapMatrix<Int> {
        val NActs = 6
        val abmMatrix = GridMapMatrix(
            IntOperators,
            (timesteps+1)*gridSize*gridSize,
            timesteps*gridSize*(gridSize-1)*NActs
        )

        // nCols = NActs*(gridsize-2)*(gridsize-2) + 2*(gridsize-2)*(NActs-2) + 2*(gridsize-2)*(NActs-1) + 4*(NActs-3)
        // =

        var actColumn = 0
        for(timestep in 0 until timesteps) {
            for(agentRow in 0 until gridSize) {
                for(agentCol in 0 until gridSize) {
                    if(agentCol>0) { // move left
                        val act = actColumn++
                        abmMatrix[eqnId(gridSize, timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(gridSize, timestep+1,agentRow,agentCol-1),act] = 1
                    }
                    if(agentCol<gridSize-1) { // move right
                        val act = actColumn++
                        abmMatrix[eqnId(gridSize, timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(gridSize,timestep+1,agentRow,agentCol+1),act] = 1
                    }
                    if(agentRow>0) { // move down
                        val act = actColumn++
                        abmMatrix[eqnId(gridSize, timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(gridSize, timestep+1,agentRow-1,agentCol),act] = 1
                    }
                    if(agentRow<gridSize-1) { // move up
                        val act = actColumn++
                        abmMatrix[eqnId(gridSize, timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(gridSize, timestep+1,agentRow+1,agentCol),act] = 1
                    }
                    if(agentCol>0) { // eat left
                        val act = actColumn++
                        abmMatrix[eqnId(gridSize, timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(gridSize, timestep,agentRow,agentCol-1),act] = -1
                        abmMatrix[eqnId(gridSize, timestep+1,agentRow,agentCol),act] = 1
                    }
                    if(agentCol<gridSize-1) { // give birth right
                        val act = actColumn++
                        abmMatrix[eqnId(gridSize, timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(gridSize, timestep+1,agentRow,agentCol+1),act] = 1
                        abmMatrix[eqnId(gridSize, timestep+1,agentRow,agentCol),act] = 1
                    }
                }
            }
        }
        return abmMatrix
    }

    fun eqnId(N:Int, timestep: Int, row: Int, col: Int) = timestep*N*N + row*N + col


//    fun twoDFermionicAbmMatrix(gridSize: Int, timesteps: Int): HashColIntMatrix {
//        val M = twoDabmMatrix(gridSize, timesteps)
//        val originalNCols = M.nCols
//        val originalNRows = M.nCols
//        M.resize(M.nRows + M.nCols,M.nCols*2)
//        for(i in originalNRows until M.nRows) {
//
//        }
//    }

    // ABM with left, right, reproduce and eat
    // T is number of timesteps
    // N is number of grid points
    fun oneDabmMatrix(T: Int, N: Int): GridMapMatrix<Int> {
        val NActs = 4
        val abmMatrix = GridMapMatrix(IntOperators, (T+1)*N, T*(N*NActs - 4))

        var actColumn = 0
        for (timestep in 0 until T) {
            for (agentPos in 0 until N) {
                if (agentPos > 0) { // move left
                    val act = actColumn++
                    abmMatrix[timestep * N + agentPos, act] = -1
                    abmMatrix[(timestep+1) * N + agentPos - 1, act] = 1
                }
                if (agentPos < N - 1) { // move right
                    val act = actColumn++
                    abmMatrix[timestep * N + agentPos, act] = -1
                    abmMatrix[(timestep+1) * N + agentPos + 1, act] = 1
                }
                if (agentPos > 0) { // eat left
                    val act = actColumn++
                    abmMatrix[timestep * N + agentPos, act] = -1
                    abmMatrix[timestep * N + agentPos - 1, act] = -1
                    abmMatrix[(timestep+1) * N + agentPos, act] = 1
                }
                if (agentPos < N - 1) { // give birth right
                    val act = actColumn++
                    abmMatrix[timestep * N + agentPos, act] = -1
                    abmMatrix[(timestep+1) * N + agentPos + 1, act] = 1
                    abmMatrix[(timestep+1) * N + agentPos, act] = 1
                }
            }
        }
        return abmMatrix
    }

    // ABM with left, right and stay-put
    // T is number of timesteps
    // N is number of grid points
    fun oneDNoInteractionMatrix(T: Int, N: Int): GridMapMatrix<Int> {
        val NActs = 3
        val abmMatrix = GridMapMatrix(IntOperators,(T+1)*N, T*(N*NActs - 2))

        var actColumn = 0
        for (timestep in 0 until T) {
            for (agentPos in 0 until N) {
                if (agentPos > 0) { // move left
                    val act = actColumn++
                    abmMatrix[timestep * N + agentPos, act] = -1
                    abmMatrix[(timestep+1) * N + agentPos - 1, act] = 1
                }
                if (agentPos < N - 1) { // move right
                    val act = actColumn++
                    abmMatrix[timestep * N + agentPos, act] = -1
                    abmMatrix[(timestep+1) * N + agentPos + 1, act] = 1
                }
                // stay put
                val act = actColumn++
                abmMatrix[timestep * N + agentPos, act] = -1
                abmMatrix[(timestep+1) * N + agentPos, act] = 1
            }
        }
        return abmMatrix
    }


    fun plotFeynmannDiagram(state: SparseIntVector, dynamics: HashColIntMatrix, agentPositions: Int, timesteps: Int) {
        val particleVectors = ArrayList<Int>() // in format (x,y,dx,dy)...
        val antiparticleVectors = ArrayList<Int>()
        for(act in state) {
            val array = if(act.value>0) particleVectors else antiparticleVectors
            val sourceLocations = ArrayList<Pair<Int,Int>>()
            val destinationLocations = ArrayList<Pair<Int,Int>>()
            for(actLocation in dynamics.columns[act.key]) {
                val location = Pair(actLocation.key.rem(agentPositions), actLocation.key.div(agentPositions))
                if(actLocation.value > 0) destinationLocations.add(location) else sourceLocations.add(location)
            }
            for(source in sourceLocations) {
                for(dest in destinationLocations) {
                    array.addAll(arrayOf(source.first, source.second, dest.first-source.first, dest.second-source.second))
                }
            }
        }

        gnuplot {
            if(particleVectors.size > 0) {
                val particleData = heredoc(particleVectors, 4)
                invoke("plot [0:${agentPositions-1}][0:$timesteps] $particleData with vectors lw 2 lc rgbcolor 0xff")
            }
            if(antiparticleVectors.size > 0) {
                val antiparticleData = heredoc(antiparticleVectors, 4)
                invoke("replot $antiparticleData with vectors lw 1 lc rgbcolor 0xff0000")
            }
        }
    }



}