package experiments

import ABMMatrices.oneDabmMatrix
import ABMMatrices.twoDabmMatrix
import HermiteDecomposition
import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import lib.SparseColIntMatrix
import lib.gnuplot
import org.apache.commons.math3.distribution.EnumeratedDistribution
import org.apache.commons.math3.primes.Primes
import org.jgrapht.alg.spanning.PrimMinimumSpanningTree
import org.jgrapht.graph.DefaultEdge
import org.jgrapht.graph.SimpleGraph
import org.junit.Test
import treeMCMC.MCMCTree
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
    fun relaxationTest() {
        // setup problem
        val abmMatrix = oneDabmMatrix(T,N)
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
    fun testStuff() {
//        val g = SparseColIntMatrix.ExtendedEuclid(1,1)
//        println("GCD=${g.gcd} x = ${g.x}  y = ${g.y}")
        val myGraph = SimpleGraph<Int,DefaultEdge>(DefaultEdge::class.java)
        myGraph.addVertex(1)
        myGraph.addVertex(2)
        myGraph.addVertex(3)
        myGraph.addVertex(4)
        myGraph.addEdge(1,2)
        myGraph.addEdge(2,3)
        myGraph.addEdge(3,4)
        myGraph.addEdge(4,1)
        myGraph.addEdge(1,3)
        myGraph.addEdge(2,4)
        println(myGraph)
        
        val spanningTree = PrimMinimumSpanningTree(myGraph).spanningTree
        println(spanningTree)
    }


    @Test
    fun MCMC() {
        // setup problem
        val abmMatrix = oneDabmMatrix(T,N)
        val b = SparseColIntMatrix.SparseIntColumn()
        b[2] = -1
        b[N*T + 3] = 1

        val mcmcTree = MCMCTree(abmMatrix, b, DoubleArray(abmMatrix.nCols) { 0.25 } )

        for(s in 1..20) {
            val sample = mcmcTree.currentSample.actValue
            println(sample)
//            plotFeynmannDiagram(sample, abmMatrix, N)
            mcmcTree.nextSample()
        }
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