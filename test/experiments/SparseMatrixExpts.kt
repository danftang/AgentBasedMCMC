package experiments

import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import hermiteForm.*
import org.apache.commons.math3.distribution.EnumeratedDistribution
import org.apache.commons.math3.primes.Primes
import org.junit.Test
import kotlin.math.ln
import kotlin.math.pow
import kotlin.random.Random

class SparseMatrixExpts {
    val N = 5
    val T = 10

    val primeCache = ArrayList<Int>()
    val lnPrimeCache = ArrayList<Double>()

    init {
        System.loadLibrary("jniortools")
        primeCache.add(2)
    }


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
        val abmMatrix = oneDabmMatrix()

        // observations
        val b = SparseColIntMatrix.SparseIntColumn()
        b[2] = -1
        b[N*T + 3] = 1

        val abm = SparseColIntMatrix(abmMatrix)
        val U = abmMatrix.hermiteDecomposition()
        val baseHermiteVector = abmMatrix.solveLowerTriangular(b)

        val baseState = U * baseHermiteVector

        // sanity check
        print("Error vector (should be empty): ")
        println(abm*baseState - b)

        val nullBasis = SparseColIntMatrix(U.subList(abmMatrix.nRows, U.nRows), U.nRows)
        nullBasis.upperTriangularise()

        val positiveDefiniteNullBasis = SparseColIntMatrix(nullBasis)
        positiveDefiniteNullBasis.upperTriangularThin()

        val shortestNullBasis = SparseColIntMatrix(nullBasis)
        shortestNullBasis.lowerTriangulariseWithoutSwap()
        for(pass in 1..10) { nullBasis.minimiseAfterDoubleTriangularisation() }

        println("Positive definite non zero elements: ${positiveDefiniteNullBasis.nonZeroElementCount()}")
        println("Shortest basis non zero elements: ${shortestNullBasis.nonZeroElementCount()}")
//        println(U)
        println(positiveDefiniteNullBasis.toSparsityString())
//        println(nullBasis)


//        nullBasis.sortBy {col ->
//            col.data.keys.max()
//        }

        val squareNullBasis = squareizeByEndPoint(positiveDefiniteNullBasis)
        val positiveBase = minimise(baseState, squareNullBasis)
        println(positiveBase)
//        plotFeynmannDiagram(positiveBase, abm, N)

//        penaltySampleShortBasis(positiveBase, shortestNullBasis)
//        penaltySamplePositiveBasis(positiveBase, positiveDefiniteNullBasis)

//        println(nullBasis.toSparsityString())
//        println(positiveBase)
//        println(abmMatrix)
//        println(U)
//        println(U[U.nCols-2])
//        println(U[U.nCols-1].toTrajectory(abm,N,N))

        // Do sampling
//        var sample = positiveBase
//        val perturbationBasis = shortestNullBasis
//        for(sampleNumber in 1..100) {
//            val possiblePerturbations = findPositiveTransitions(sample, perturbationBasis)
//            println("No of perturbations ${possiblePerturbations.size}")
//            if(possiblePerturbations.isEmpty()) {
//                println("No transitions")
//            } else {
//                sample = possiblePerturbations[Random.nextInt(possiblePerturbations.size)]
//                println(sample)
//            }
//        }


//        println(baseState)
        plotFeynmannDiagram(positiveDefiniteNullBasis[67], abm, N)
        plotFeynmannDiagram(positiveDefiniteNullBasis[68], abm, N)
        plotFeynmannDiagram(positiveDefiniteNullBasis[69], abm, N)
        plotFeynmannDiagram(positiveDefiniteNullBasis[70], abm, N)
//        plotFeynmannDiagram(nullBasis[1], abm, N)
//        plotFeynmannDiagram(nullBasis[2], abm, N)
//        plotFeynmannDiagram(nullBasis[nullBasis.nCols-1], abm, N)
//        plotFeynmannDiagram(nullBasis[nullBasis.nCols-2], abm, N)
//        plotFeynmannDiagram(nullBasis[nullBasis.nCols-3], abm, N)
//        plotTrajectory(possiblePerturbations[1], abm,N,N)

//        printAllActs(abm,N,N)
    }

    @Test
    fun relaxationTest() {
        // setup problem
        val abmMatrix = oneDabmMatrix()
        val b = SparseColIntMatrix.SparseIntColumn()
        b[2] = -1
        b[N*T + 3] = 1

        var prime = 1
        val nBases = abmMatrix.nCols
        val firstNPrimes = IntArray(nBases) {i ->
            prime = Primes.nextPrime(prime+1)
            prime
        }
        val objectiveCoeffs = Array(nBases) {i ->
            ln(firstNPrimes[nBases - i - 1].toDouble())
        }

        val solution = LPsolve(abmMatrix, b, objectiveCoeffs)
        println(solution.asList())
    }

    @Test
    fun findSmallestSolution() {

        // setup problem
        val abmMatrix = oneDabmMatrix()
        val abm = SparseColIntMatrix(abmMatrix)
        val b = SparseColIntMatrix.SparseIntColumn()
        b[2] = -1
        b[N*T + 3] = 1

        // Hermite decompose
        val U = abmMatrix.hermiteDecomposition()
        val baseHermiteVector = abmMatrix.solveLowerTriangular(b)
        val baseState = U * baseHermiteVector
        val nullBasis = SparseColIntMatrix(U.subList(abmMatrix.nRows, U.nRows), U.nRows)
        nullBasis.upperTriangularise()

        val positiveDefiniteNullBasis = SparseColIntMatrix(nullBasis)
        positiveDefiniteNullBasis.upperTriangularThin()
        val positiveBaseState = upperTriangularMinimise(baseState, positiveDefiniteNullBasis)
        println(positiveDefiniteNullBasis.toSparsityString())
        println(generateBasisConstraintMatrix(positiveDefiniteNullBasis).transpose().toSparsityString())
        println(positiveBaseState)

//        val smallestVec = IPsolve(positiveDefiniteNullBasis, -positiveBaseState, getReverseLogPrimeVector(positiveDefiniteNullBasis.nCols))
//        val basisConstraints = generateBasisConstraintMatrix(positiveDefiniteNullBasis).transpose()
//        val smallestVec = IPsolve(basisConstraints, -positiveBaseState, DoubleArray(positiveDefiniteNullBasis.nCols) {0.0}.asList())
        val smallestVec = IPsolve(positiveDefiniteNullBasis, -positiveBaseState, DoubleArray(positiveDefiniteNullBasis.nCols) {0.0}.asList())
        println(smallestVec.asList())
        val solution = positiveBaseState + positiveDefiniteNullBasis * SparseColIntMatrix.SparseIntColumn(smallestVec)
        plotFeynmannDiagram(solution,abm, N)

    }

    @Test
    fun MCMC() {
        // setup problem
        val abmMatrix = oneDabmMatrix()
        val b = SparseColIntMatrix.SparseIntColumn()
        b[2] = -1
        b[N*T + 3] = 1

        val mcmcTree = MCMCTree(abmMatrix, b, DoubleArray(abmMatrix.nCols) { 0.25 } )

        for(s in 1..3) {
            val sample = mcmcTree.currentSample.actValue
            println(sample)
            plotFeynmannDiagram(sample, abmMatrix, N)
            mcmcTree.nextSample()
        }
//        println(abmMatrix.toSparsityString())

//        val decomposition = HermiteDecomposition(abmMatrix, b)
//        println(decomposition.nullBasis.toSparsityString())
//        val baseSolution = decomposition.baseSolution(b)

//        val constraints = decomposition.toConstraintMatrix()
//        val constraintVector = decomposition.actToNonBasisVector(-decomposition.baseState)
//        val rootBasisVector = constraints.IPsolve(constraintVector, DoubleArray(constraints.nCols) {0.0}.asList())
//        val rootSolution = decomposition.basisVectorToActs(rootBasisVector)
//        plotFeynmannDiagram(rootSolution, abmMatrix, N)

//        println(rootBasisVector.asList())
//        println(rootBasisVector[34])


//        var sample = MCMCTreeNode(rootBasisVector, rootBasisVector.lastIndex, null)

//        println("Root sample is ${decomposition.basisVectorToActs(sample.value)}")

//        sample.generateAllChildren(decomposition)
//        sample.children.forEach {
//            val solution = decomposition.basisVectorToActs(it.value)
////            plotFeynmannDiagram(solution, abmMatrix, N)
//        }
//        for(child in sample.children) {
//            child.generateAllChildren(decomposition)
//        }



//            for (n in 1 until sample.breakpoint) {
//                try {
//                    val newBreakpoint = n
//                    if (newBreakpoint != null) {
//                        println("Trying Break at $newBreakpoint")
//                        val childSample = sample.splitAt(newBreakpoint, decomposition)
//                        val solution = decomposition.basisVectorToActs(childSample.value)
//                        println("Sample: $solution")
////                    plotFeynmannDiagram(solution, abmMatrix, N)
//                    }
//                } catch (err: java.lang.RuntimeException) {
////                    println("Got runtime exception")
//                }
//            }

    }

    fun generateBasisConstraintMatrix(positiveDefiniteNullBasis: SparseColIntMatrix): SparseColIntMatrix {
        val constraints = positiveDefiniteNullBasis.transpose()
        constraints.removeIf { it.entries.size <= 1 }
        return constraints
    }

    // true if the act with this index is in the basis
    fun getPositiveDefiniteBasisMask(positiveDefiniteNullBasis: SparseColIntMatrix): BooleanArray {
        val mask = BooleanArray(positiveDefiniteNullBasis.nRows) { false }
        for(col in positiveDefiniteNullBasis.columns) {
            val basisAct = col.keys.max()
            if(basisAct != null) mask[basisAct] = true
        }
        return mask
    }

    fun getReverseLogPrimeVector(size: Int): List<Double> {
        while(primeCache.size < size) {
            val nextPrime = Primes.nextPrime(primeCache.last()+1)
            primeCache.add(nextPrime)
            lnPrimeCache.add(ln(nextPrime.toDouble()))
        }
        return lnPrimeCache.subList(0, size).asReversed()
    }

    // Minimise CX
    // Subject to:
    // MX >= B
    // and
    // X >= 0 and X is integer
    //
    fun IPsolve(M: SparseColIntMatrix, B: SparseColIntMatrix.SparseIntColumn, C: List<Double>): IntArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val X = solver.makeIntVarArray(M.nCols, 0.0, Double.POSITIVE_INFINITY)
        val constraints = Array<MPConstraint>(M.nRows) { solver.makeConstraint() }
        for(col in 0 until M.nCols) {
            for(coefficient in M[col].entries) {
                constraints[coefficient.key].setCoefficient(X[col], coefficient.value.toDouble())
            }
        }
        for(i in 0 until M.nRows) {
            constraints[i].setBounds(B[i].toDouble(), Double.POSITIVE_INFINITY)
        }
        val objective = solver.objective()
        for(i in C.indices) {
            objective.setCoefficient(X[i], C[i])
        }
        solver.objective().setMinimization()
        val solveState = solver.solve()
        return if (solveState == MPSolver.ResultStatus.OPTIMAL)
            IntArray(X.size) { i -> X[i].solutionValue().toInt() }
        else
            throw(RuntimeException(
                when (solveState) {
                    MPSolver.ResultStatus.INFEASIBLE -> "Infeasible"
                    MPSolver.ResultStatus.UNBOUNDED -> "Unbounded"
                    else -> "Solve Error"
                }
            ))
    }


    // Minimise CX
    // Subject to:
    // MX >= B
    // and
    // X >= 0
    //
    fun LPsolve(M: SparseColIntMatrix, B: SparseColIntMatrix.SparseIntColumn, C: Array<Double>): DoubleArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.GLOP_LINEAR_PROGRAMMING)
        val X = solver.makeNumVarArray(M.nCols, 0.0, Double.POSITIVE_INFINITY)
        val constraints = Array<MPConstraint>(M.nRows) { solver.makeConstraint() }
        for(col in 0 until M.nCols) {
            for(coefficient in M[col].entries) {
                constraints[coefficient.key].setCoefficient(X[col], coefficient.value.toDouble())
            }
        }
        for(i in 0 until M.nRows) {
            constraints[i].setBounds(B[i].toDouble(), Double.POSITIVE_INFINITY)
        }
        val objective = solver.objective()
        for(i in C.indices) {
            objective.setCoefficient(X[i], C[i])
        }
        solver.objective().setMinimization()
        val solveState = solver.solve()
        return if (solveState == MPSolver.ResultStatus.OPTIMAL)
            DoubleArray(X.size) { i -> X[i].solutionValue() }
        else
            throw(RuntimeException(
                when (solveState) {
                    MPSolver.ResultStatus.INFEASIBLE -> "Infeasible"
                    MPSolver.ResultStatus.UNBOUNDED -> "Unbounded"
                    else -> "Solve Error"
                }
            ))
    }


    // Minimise CX
    // Subject to:
    // MX >= B
    // and
    // X >= 0
    // and
    // Xi is an integer if it is an interaction
    //
    fun MIPsolve(M: SparseColIntMatrix, B: SparseColIntMatrix.SparseIntColumn, C: Array<Double>): DoubleArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val X = Array<MPVariable>(M.nCols) {i ->
            if(M[i].entries.size > 2) {
                // interaction
                solver.makeIntVar(0.0,Double.POSITIVE_INFINITY,"I$i")
            } else {
                println("Making double var")
                solver.makeNumVar(0.0,Double.POSITIVE_INFINITY,"D$i")
            }
        }
        val constraints = Array<MPConstraint>(M.nRows) { solver.makeConstraint() }
        for(col in 0 until M.nCols) {
            for(coefficient in M[col].entries) {
                constraints[coefficient.key].setCoefficient(X[col], coefficient.value.toDouble())
            }
        }
        for(i in 0 until M.nRows) {
            constraints[i].setBounds(B[i].toDouble(), Double.POSITIVE_INFINITY)
        }
        val objective = solver.objective()
        for(i in C.indices) {
            objective.setCoefficient(X[i], C[i])
        }
        solver.objective().setMinimization()
        val solveState = solver.solve()
        return if (solveState == MPSolver.ResultStatus.OPTIMAL)
            DoubleArray(X.size) { i -> X[i].solutionValue() }
        else
            throw(RuntimeException(
                when (solveState) {
                    MPSolver.ResultStatus.INFEASIBLE -> "Infeasible"
                    MPSolver.ResultStatus.UNBOUNDED -> "Unbounded"
                    else -> "Solve Error"
                }
            ))
    }


    fun findPositiveTransitions(startState: SparseColIntMatrix.SparseIntColumn, nullBasis: SparseColIntMatrix): List<SparseColIntMatrix.SparseIntColumn> {
        val transitions = ArrayList<SparseColIntMatrix.SparseIntColumn>()
        for(j in 0 until nullBasis.nCols) {
            val suggestion = startState + nullBasis[j]
            if(suggestion.isPositive()) {
                transitions.add(suggestion)
            } else {
                val suggestion2 = startState - nullBasis[j]
                if (suggestion2.isPositive()) transitions.add(suggestion2)
            }
        }
        return transitions
    }

    fun penaltySampleShortBasis(initialSample: SparseColIntMatrix.SparseIntColumn, shortNullBasis: SparseColIntMatrix) {
        val perturbationBasis = SparseColIntMatrix(shortNullBasis.nRows, shortNullBasis.nCols*2)
        for(i in 0 until shortNullBasis.nCols) {
            perturbationBasis[2*i] = shortNullBasis[i]
            perturbationBasis[2*i+1] = -shortNullBasis[i]
        }
        var currentSample = initialSample
        for(sample in 1..1000) {
            currentSample = perturbWithPenalty(currentSample, perturbationBasis)
            if(currentSample.isPositive()) println(currentSample) else println(1.0/negativePenalty(currentSample))
        }
    }

    fun penaltySamplePositiveBasis(initialSample: SparseColIntMatrix.SparseIntColumn, positiveNullBasis: SparseColIntMatrix) {
        val nullToActMap = Array(positiveNullBasis.nCols) { col ->
            positiveNullBasis[col].keys.max()?:0
        }
        var currentSample = initialSample
        for(sample in 1..10000) {
            currentSample = perturbWithPenalty(currentSample, positiveNullBasis, nullToActMap)
            if(currentSample.isPositive()) println(currentSample) // else println(1.0/negativePenalty(currentSample))
        }
    }


    fun perturbWithPenalty(currentSample: SparseColIntMatrix.SparseIntColumn, positiveNullBasis: SparseColIntMatrix, nullToActBasisIndexMap: Array<Int>): SparseColIntMatrix.SparseIntColumn {
        val pmf = ArrayList<org.apache.commons.math3.util.Pair<SparseColIntMatrix.SparseIntColumn, Double>>(positiveNullBasis.nCols)
        for(basisIndex in 0 until positiveNullBasis.nCols) {
            val currentbasisVal = currentSample[nullToActBasisIndexMap[basisIndex]]
            val newSample = if(currentbasisVal == 0) {
                currentSample + positiveNullBasis[basisIndex]
            } else {
                currentSample - positiveNullBasis[basisIndex]
            }
            pmf.add(org.apache.commons.math3.util.Pair(newSample, negativePenalty(newSample)))
        }
        val choose = EnumeratedDistribution(pmf)
        return choose.sample()
    }

    fun perturbWithPenalty(currentSample: SparseColIntMatrix.SparseIntColumn, perturbations: SparseColIntMatrix): SparseColIntMatrix.SparseIntColumn {
        val pmf = ArrayList<org.apache.commons.math3.util.Pair<SparseColIntMatrix.SparseIntColumn, Double>>(perturbations.nCols)
        for(proposal in perturbations) {
            val newSample = currentSample + proposal
            pmf.add(org.apache.commons.math3.util.Pair(newSample, negativePenalty(newSample)))
        }
        val choose = EnumeratedDistribution(pmf)
        return choose.sample()
    }

    fun negativePenalty(basisVector: SparseColIntMatrix.SparseIntColumn): Double {
        val penaltyPerNegative = 0.5
        val negativeCount = basisVector.values.sumBy { if(it < 0) -it else 0 }
        return penaltyPerNegative.pow(negativeCount)
    }



    fun twoDabmMatrix(): SparseColIntMatrix {
        val NActs = 6
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

    fun oneDabmMatrix(): SparseColIntMatrix {
        val NActs = 4
        val abmMatrix = SparseColIntMatrix((T+1)*N, T*(N*NActs - 4))

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


    fun eqnId(timestep: Int, row: Int, col: Int) = timestep*N*N + row*N + col

    fun squareizeByEndPoint(matrix: SparseColIntMatrix): SparseColIntMatrix {
        val result = SparseColIntMatrix(matrix.nRows, matrix.nRows)
        for(column in matrix.columns) {
            val endPoint = column.data.keys.max()?:0
            result.columns[endPoint] = column
        }
        return result
    }

    fun minimise(col: SparseColIntMatrix.SparseIntColumn, squareNullBasis: SparseColIntMatrix): SparseColIntMatrix.SparseIntColumn {
        val denseCol = col.toIntArray(squareNullBasis.nCols)
        for(i in denseCol.size-1 downTo 0) {
            if(denseCol[i] != 0 && squareNullBasis[i,i] != 0) {
                val w = -denseCol[i] / squareNullBasis[i, i]
                for (element in squareNullBasis[i].entries) {
                    denseCol[element.key] += w * element.value
                }
            }
        }
        return SparseColIntMatrix.SparseIntColumn(denseCol)
    }

    fun upperTriangularMinimise(toMinimise: SparseColIntMatrix.SparseIntColumn, upperNullBasis: SparseColIntMatrix): SparseColIntMatrix.SparseIntColumn {
        val minimised = toMinimise.copy()
        for(basis in upperNullBasis.asReversed()) {
            val maxIndex = basis.keys.max()
            if(maxIndex != null) {
                minimised.weightedPlusAssign(basis, -minimised[maxIndex]/basis[maxIndex])
            }
        }
        return minimised
    }

    fun plot2DTrajectory(column: SparseColIntMatrix.SparseIntColumn, dynamics: SparseColIntMatrix, agentRows: Int, agentCols: Int) {
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
            if(particleVectors.size > 0) {
                val particleData = heredoc(particleVectors, 6)
                invoke("splot [0:$T][0:$N][0:$N] $particleData with vectors lw 2 lc rgbcolor 0xff")
            }
            if(antiparticleVectors.size > 0) {
                val antiparticleData = heredoc(antiparticleVectors, 6)
                invoke("replot $antiparticleData with vectors lw 1 lc rgbcolor 0xff0000")
            }
        }
    }


    fun plotFeynmannDiagram(column: SparseColIntMatrix.SparseIntColumn, dynamics: SparseColIntMatrix, agentPositions: Int) {
        val particleVectors = ArrayList<Int>() // in format (x,y,dx,dy)...
        val antiparticleVectors = ArrayList<Int>()
        for(act in column.entries) {
            val array = if(act.value>0) particleVectors else antiparticleVectors
            val sourceLocations = ArrayList<Pair<Int,Int>>()
            val destinationLocations = ArrayList<Pair<Int,Int>>()
            for(actLocation in dynamics[act.key].entries) {
                val location = Pair(actLocation.key.rem(N), actLocation.key.div(N))
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
                invoke("plot [0:${N-1}][0:$T] $particleData with vectors lw 2 lc rgbcolor 0xff")
            }
            if(antiparticleVectors.size > 0) {
                val antiparticleData = heredoc(antiparticleVectors, 4)
                invoke("replot $antiparticleData with vectors lw 1 lc rgbcolor 0xff0000")
            }
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