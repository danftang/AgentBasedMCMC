package experiments

import ABMMatrices.oneDabmMatrix
import ABMMatrices.plotFeynmannDiagram
import ABMMatrices.twoDabmMatrix
import ConvexPolyhedron
import lib.sparseMatrix.HashColIntMatrix
import lib.sparseMatrix.HashIntVector
import lib.sparseMatrix.HashRowColIntMatrix
import lib.sparseMatrix.SparseIntVector
import org.junit.Test
import kotlin.random.Random

class HermiteFormExpts {

    init {
        System.loadLibrary("jniortools")
    }

    @Test
    fun testABMtoHermite() {
        val T = 3
        val N = 5
        val abmMatrix = oneDabmMatrix(T,N)

        // observations
        val b = HashIntVector()
        b[2] = -1
        b[N*T + 3] = 1

        val abm = HashColIntMatrix(abmMatrix)
        val U = abmMatrix.hermiteDecomposition()
        val baseHermiteVector = abmMatrix.solveLowerTriangular(b)

        val baseState = U * baseHermiteVector

        // sanity check
        print("Error vector (should be empty): ")
        println(abm*baseState - b)

        val nullBasis = HashColIntMatrix(U.subMatrixView(abmMatrix.nRows, U.nRows))
        nullBasis.upperTriangularise()

        val positiveDefiniteNullBasis = HashColIntMatrix(nullBasis)
        positiveDefiniteNullBasis.upperTriangularThin()

        val shortestNullBasis = HashColIntMatrix(nullBasis)
        shortestNullBasis.lowerTriangulariseWithoutSwap()
        for(pass in 1..10) { nullBasis.minimiseAfterDoubleTriangularisation() }

        println("Positive definite non zero elements: ${positiveDefiniteNullBasis.nonZeroElementCount()}")
        println("Shortest basis non zero elements: ${shortestNullBasis.nonZeroElementCount()}")
//        println(U)
        println(positiveDefiniteNullBasis.toSparsityString())
//        println(nullBasis)

        val squareNullBasis = squareizeByEndPoint(positiveDefiniteNullBasis)
        val positiveBase = minimise(baseState, squareNullBasis)
        println(positiveBase)
        plotFeynmannDiagram(positiveDefiniteNullBasis.columns[67], abm, N, T)
        plotFeynmannDiagram(positiveDefiniteNullBasis.columns[68], abm, N, T)
    }

    @Test
    fun treeDecomposition() {
        val gridSize = 32
        val timesteps = 6
        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
        println("abm sparsity is ${abmMatrix.sparsityRatio()}")
        println("abm geometry is ${abmMatrix.nRows} x ${abmMatrix.nCols}")
//        println(abmMatrix.toSparsityString())
        val observations = HashIntVector()
        for(agent in 1..20) {
            val xPos = Random.nextInt(gridSize)
            val yPos = Random.nextInt(gridSize)
            observations[xPos*gridSize + yPos] = -1
            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
        }

        // Solve ABMMatrix directly
        println("Starting direct solve")
//        val max = solution.entries.first()
//        println("trying to maximise ${max.key}")
//        val positiveDirectSolution = abmMatrix.IPsolve(observations, DoubleArray(abmMatrix.nCols) { if(it == max.key) -1.0 else 0.0 }.asList(),"=")
        val positiveDirectSolution = abmMatrix.IPsolve(observations, DoubleArray(abmMatrix.nCols) {0.0}.asList(), "=")
        println(positiveDirectSolution)
        println("Error vector (should be empty) = ${abmMatrix * positiveDirectSolution - observations}")


        // zero removal
        val polyhedron = ConvexPolyhedron(HashRowColIntMatrix(abmMatrix), HashIntVector(observations))
        println("Original polyhedron size is ${polyhedron.M.nRows} ${polyhedron.M.nCols}")
        polyhedron.removeZeros()
        println("Reduced polyhedron size is ${polyhedron.M.nRows} ${polyhedron.M.nCols}")

        // direct solve after zero removal
        println("Starting reduced direct solve")
        val positiveDirectSolution2 = polyhedron.findValidPoint()
        println(positiveDirectSolution2)
        println("Error vector (should be empty) = ${abmMatrix * positiveDirectSolution2 - observations}")


        // Hermite decompose
//        val rootRow =  polyhedron.M.nRows / 2 //gridSize*(gridSize*(timesteps/2) + gridSize/2) + gridSize/2
//        println("Creating tree at root row $rootRow")
//        val hermite = HermiteDecomposition(polyhedron.M, polyhedron.B, rootRow)
//        println("H sparsity is ${hermite.H.sparsityRatio()}")
//        println("H range = ${hermite.H.valueRange()}")
////        println(hermite.H.toSparsityString())
////        println(hermite.H)
//        println("U sparsity is ${hermite.U.sparsityRatio()}")
//        println("U range = ${hermite.U.valueRange()}")
////        println(hermite.U.toSparsityString())
////        println(hermite.U)
//
//        assert(hermite.nullBasis.checkNoZeroEntries())
//        println("Solution basis range is ${hermite.solutionBasis.valueRange()}")
//        println("null value range is ${hermite.nullBasis.valueRange()}")
//        println("Null basis sparsity is ${hermite.nullBasis.sparsityRatio()}")
//
////        println("abmMatrix size is ${abmMatrix.nRows} ${abmMatrix.nCols}")
//
//        println("Error vector (should be empty) = ${abmMatrix * T * hermite.baseState - observations}")
//
//        println("Starting null basis solve")
//        val reducedNullBasis = hermite.nullMaskMatrix * hermite.nullBasis
//        val reducedBaseState = hermite.nullMaskMatrix * hermite.baseState
//        val positiveBasisSolution = reducedNullBasis.IPsolve(-reducedBaseState, DoubleArray(hermite.nullBasis.nCols) {0.0}.asList())
//        val solution = hermite.baseState + hermite.nullBasis * positiveBasisSolution
//        println(solution)


    }

    @Test
    fun scaleUp() {
        val gridSize = 6
        val timesteps = 2
        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
        println("abm sparsity is ${abmMatrix.sparsityRatio()}")
//        println(abmMatrix.toSparsityString())
        val observations = HashIntVector()
        for(agent in 1..20) {
            val xPos = Random.nextInt(gridSize)
            val yPos = Random.nextInt(gridSize)
            observations[xPos*gridSize + yPos] = -1
            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
        }

        // Hermite decompose
        val hermiteForm = abmMatrix.copy()
        println("Starting...")
        val U = hermiteForm.hermiteDecomposition()
//        println(hermiteForm.diagonal())
//        println(U.toSparsityString())
//        println(hermiteForm.toSparsityString())
        println("Done Hermite decomposition")
        // AX = B, X = UY, HY = B
        val baseHermiteVector = hermiteForm.solveLowerTriangular(observations)
        println("Base hermite state is $baseHermiteVector")
        val baseState = U * baseHermiteVector
        println("Base state is $baseState")
        println("Hermite Residual is ${hermiteForm * baseHermiteVector - observations}")
        println("Act Residual is ${abmMatrix * baseState - observations}")
        val nullBasis = U.subMatrixView(hermiteForm.nRows, U.nCols)
        val solutionBasis = U.subMatrixView(0, hermiteForm.nRows)
        println("Solved for base state and generated null basis")
//        println(nullBasis)

        assert(nullBasis.checkNoZeroEntries())
        println("Solution basis range is ${solutionBasis.valueRange()}")
        println("null value range is ${nullBasis.valueRange()}")
        println("Null basis sparsity is ${nullBasis.sparsityRatio()}")

        nullBasis.upperTriangularise()
        assert(nullBasis.checkNoZeroEntries())
        println("Upper triangularised null basis")
        println("null value range is ${nullBasis.valueRange()}")
        println("Null basis sparsity is ${nullBasis.sparsityRatio()}")
//        println(nullBasis)
        nullBasis.upperTriangularThin()
        assert(nullBasis.checkNoZeroEntries())
        println("Thinned null basis")
        println("null value range is ${nullBasis.valueRange()}")
        println("Null basis sparsity is ${nullBasis.sparsityRatio()}")
        val positiveDefiniteBaseState = upperTriangularMinimise(baseState, nullBasis)
        println("Got positive definite null basis and base state")
        println(positiveDefiniteBaseState)


        val positiveBasisSolution = nullBasis.IPsolve(-positiveDefiniteBaseState, DoubleArray(nullBasis.nCols) {0.0}.asList())
        println("Solved for positive solution")
        val solution = positiveDefiniteBaseState + nullBasis * positiveBasisSolution
        println(solution)
    }

    @Test
    fun findSmallestSolution() {

        // setup problem
        val T = 2
        val N = 5
        val abmMatrix = oneDabmMatrix(T,N)
        val abm = HashColIntMatrix(abmMatrix)
        val b = HashIntVector()
        b[2] = -1
        b[N*T + 3] = 1

        // Hermite decompose
        val U = abmMatrix.hermiteDecomposition()
        val baseHermiteVector = abmMatrix.solveLowerTriangular(b)
        val baseState = U * baseHermiteVector
        val nullBasis = HashColIntMatrix(U.subMatrixView(abmMatrix.nRows, U.nRows))
        nullBasis.upperTriangularise()

        val positiveDefiniteNullBasis = HashColIntMatrix(nullBasis)
        positiveDefiniteNullBasis.upperTriangularThin()
        val positiveBaseState = upperTriangularMinimise(baseState, positiveDefiniteNullBasis)
        println(positiveDefiniteNullBasis.toSparsityString())
//        println(generateBasisConstraintMatrix(positiveDefiniteNullBasis).transpose().toSparsityString())
        println(positiveBaseState)

//        val smallestVec = IPsolve(positiveDefiniteNullBasis, -positiveBaseState, getReverseLogPrimeVector(positiveDefiniteNullBasis.nCols))
//        val basisConstraints = generateBasisConstraintMatrix(positiveDefiniteNullBasis).transpose()
//        val smallestVec = IPsolve(basisConstraints, -positiveBaseState, DoubleArray(positiveDefiniteNullBasis.nCols) {0.0}.asList())
        val smallestVec = positiveDefiniteNullBasis.IPsolve(-positiveBaseState, DoubleArray(positiveDefiniteNullBasis.nCols) {0.0}.asList())
        println(smallestVec)
        val solution = positiveBaseState + positiveDefiniteNullBasis * smallestVec
        plotFeynmannDiagram(solution.toSparseIntVector(),abm, N, T)
    }


    fun upperTriangularMinimise(toMinimise: SparseIntVector, upperNullBasis: HashColIntMatrix): HashIntVector {
        val minimised = HashIntVector(toMinimise)
        for(basis in upperNullBasis.columns.asReversed()) {
            val maxIndex = basis.keys.max()
            if(maxIndex != null) {
                if(basis[maxIndex] == 0) println("Found zero entry in column $basis")
                minimised.weightedPlusAssign(basis, -minimised[maxIndex]/basis[maxIndex])
            }
        }
        return minimised
    }

    fun squareizeByEndPoint(matrix: HashColIntMatrix): HashColIntMatrix {
        val result = HashColIntMatrix(matrix.nRows, matrix.nRows)
        for(column in matrix.columns) {
            val endPoint = column.keys.max()?:0
            result.setColumn(endPoint, column)
        }
        return result
    }

    fun minimise(col: SparseIntVector, squareNullBasis: HashColIntMatrix): SparseIntVector {
        val denseCol = col.toIntVector(squareNullBasis.nCols)
        for(i in denseCol.size-1 downTo 0) {
            if(denseCol[i] != 0 && squareNullBasis[i,i] != 0) {
                val w = -denseCol[i] / squareNullBasis[i, i]
                for (element in squareNullBasis.columns[i]) {
                    denseCol[element.key] += w * element.value
                }
            }
        }
        return denseCol.toSparseIntVector()
    }



}