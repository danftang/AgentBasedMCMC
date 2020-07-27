package experiments

import cnf.*
import hermiteForm.SparseColIntMatrix
import hermiteForm.gnuplot
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D
import org.junit.Test
import org.ujmp.core.SparseMatrix
import org.ujmp.core.intmatrix.SparseIntMatrix
import org.ujmp.core.intmatrix.SparseIntMatrix2D
import org.ujmp.core.matrix.factory.SparseMatrixFactory

class SparseMatrixExpts {
    val N = 3
    val T = 3
    val NActs = 6


    @Test
    fun simpleOpsTest() {
        val mySparseMatrix = SparseColIntMatrix.Identity(10)
        mySparseMatrix[3,6] = 4
        mySparseMatrix[0].weightedPlusAssign(mySparseMatrix[6],-2)
        mySparseMatrix[6].weightedPlusAssign(mySparseMatrix[0],1)
        println(mySparseMatrix)
    }


    @Test
    fun testABMtoHermite() {
//        val abmMatrix = SparseColIntMatrix((T+1)*N*N, T*N*N*NActs)
//
//        var lastTimestepEqns = Array(N) { Array(N) { HashMap<Int,Int>() } }
//        var nextTimestepEqns = Array(N) { Array(N) { HashMap<Int,Int>() } }
//
//        var actColumn = 0
//        var eqnRow = 0
//        for(t in 1..T) {
//            for (agentRow in 0 until N) {
//                for (agentCol in 0 until N) {
//                    val left = actColumn++
//                    val right = actColumn++
//                    val up = actColumn++
//                    val down = actColumn++
//                    val eatLeft = actColumn++
//                    val giveBirth = actColumn++
//
//                    val eqn = lastTimestepEqns[agentRow][agentCol]
//                    eqn[left] = -1
//                    eqn[right] = -1
//                    eqn[up] = -1
//                    eqn[down] = -1
//                    eqn[giveBirth] = -1
//                    eqn[eatLeft] = -1
//                    lastTimestepEqns[agentRow][(agentCol + N - 1).rem(N)][eatLeft] = -1
//
//
//                    nextTimestepEqns[agentRow][(agentCol + N - 1).rem(N)][left] = 1
//                    nextTimestepEqns[agentRow][(agentCol + 1).rem(N)][right] = 1
//                    nextTimestepEqns[(agentRow + 1).rem(N)][agentCol][up] = 1
//                    nextTimestepEqns[(agentRow + N - 1).rem(N)][agentCol][down] = 1
//                    nextTimestepEqns[agentRow][agentCol][eatLeft] = 1
//                    nextTimestepEqns[agentRow][agentCol][giveBirth] = 1
//                    nextTimestepEqns[agentRow][(agentCol+1).rem(N)][giveBirth] = 1
//
//                    abmMatrix.setRow(eqnRow++, lastTimestepEqns[agentRow][agentCol].entries)
//                }
//            }
//            lastTimestepEqns = nextTimestepEqns
//            nextTimestepEqns = Array(N) { Array(N) { HashMap<Int,Int>() } }
//        }
//        for (agentRow in 0 until N) {
//            for (agentCol in 0 until N) {
//                abmMatrix.setRow(eqnRow++, lastTimestepEqns[agentRow][agentCol].entries)
//            }
//        }

        val abmMatrix = abmMatrix()

        // observations
        val b = HashMap<Int,Int>()
        var startStateRow = 0
        var endStateRow = N*N*(T-1)
        for (row in 0 until N) {
            for (col in 0 until N) {
                b[startStateRow++] = if(row == 1 && col == 1)  1 else 0
                b[endStateRow++] = if(row == 0 && col == 0) 1 else 0
            }
        }

        println(abmMatrix.toSparsityString())

        val abm = SparseColIntMatrix(abmMatrix)
        val U = abmMatrix.hermiteDecomposition()

        println(abmMatrix.toSparsityString())
        println(U.toSparsityString())

//        println(abmMatrix)
//        println(U)

//        println(U[U.nCols-2])
        println(U[U.nCols-1].toTrajectory(abm,N,N))
        plotTrajectory(U[U.nCols-1], abm,N,N)

        printAllActs(abm,N,N)
    }


    fun abmMatrix(): SparseColIntMatrix {
        val abmMatrix = SparseColIntMatrix((T+1)*N*N, T*N*(N-1)*NActs)

        var actColumn = 0
        for(timestep in 0 until T) {
            for(agentRow in 0 until N) {
                for(agentCol in 0 until N) {
                    if(agentCol>0) { // move left
                        val act = actColumn++
                        abmMatrix[eqnId(timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(timestep+1,agentRow,agentCol-1),act] = 1
                    }
                    if(agentCol<N-1) { // move right
                        val act = actColumn++
                        abmMatrix[eqnId(timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(timestep+1,agentRow,agentCol+1),act] = 1
                    }
                    if(agentRow>0) { // move down
                        val act = actColumn++
                        abmMatrix[eqnId(timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(timestep+1,agentRow-1,agentCol),act] = 1
                    }
                    if(agentRow<N-1) { // move up
                        val act = actColumn++
                        abmMatrix[eqnId(timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(timestep+1,agentRow+1,agentCol),act] = 1
                    }
                    if(agentCol>0) { // eat left
                        val act = actColumn++
                        abmMatrix[eqnId(timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(timestep,agentRow,agentCol-1),act] = -1
                        abmMatrix[eqnId(timestep+1,agentRow,agentCol),act] = 1
                    }
                    if(agentCol<N-1) { // give birth right
                        val act = actColumn++
                        abmMatrix[eqnId(timestep,agentRow,agentCol),act] = -1
                        abmMatrix[eqnId(timestep+1,agentRow,agentCol+1),act] = 1
                        abmMatrix[eqnId(timestep+1,agentRow,agentCol),act] = 1
                    }
                }
            }
        }
        return abmMatrix
    }


    fun eqnId(timestep: Int, row: Int, col: Int) = timestep*N*N + row*N + col


    fun plotTrajectory(column: SparseColIntMatrix.SparseIntColumn, dynamics: SparseColIntMatrix, agentRows: Int, agentCols: Int) {
        val particleVectors = ArrayList<Int>()
        val antiparticleVectors = ArrayList<Int>()
        for(act in column.entries) {
            val array = if(act.value>0) particleVectors else antiparticleVectors
            val sourceLocations = ArrayList<Triple<Int,Int,Int>>()
            val destinationLocations = ArrayList<Triple<Int,Int,Int>>()
            for(actLocation in dynamics[act.key].entries) {
                val location = Triple(actLocation.key.div(N*N), actLocation.key.rem(agentCols), actLocation.key.div(agentCols).rem(agentRows))
                if(actLocation.value > 0) destinationLocations.add(location) else sourceLocations.add(location)
            }
            for(source in sourceLocations) {
                for(dest in destinationLocations) {
                    array.addAll(arrayOf(source.first, source.second, source.third, dest.first-source.first, dest.second-source.second, dest.third-source.third))
                }
            }
            println()
        }

        gnuplot {
            val particleData = heredoc(particleVectors, 6)
            val antiparticleData = heredoc(antiparticleVectors, 6)
            invoke("splot $particleData with vectors lw 2 lc rgbcolor 0xff")
            invoke("replot $antiparticleData with vectors lw 2 lc rgbcolor 0xff0000")
        }
    }


    fun printAllActs(dynamics: SparseColIntMatrix, agentCols: Int, agentRows: Int) {
        for (act in dynamics) {
            val sourceLocations = ArrayList<Pair<Int, Int>>()
            val destinationLocations = ArrayList<Pair<Int, Int>>()
            for (actLocation in act.entries) {
                val location = Pair(actLocation.key.rem(agentCols), actLocation.key.div(agentCols).rem(agentRows))
                if (actLocation.value > 0) destinationLocations.add(location) else sourceLocations.add(location)
            }
            for (source in sourceLocations) {
                for (dest in destinationLocations) {
                    println("(${source.first},${source.second})->(${dest.first},${dest.second})")
                }
            }
            println()
        }
    }

}