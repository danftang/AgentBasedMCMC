import com.google.ortools.linearsolver.MPSolver

import ORTools.setConstraintsWithExplicitSlacks
import com.google.ortools.linearsolver.MPObjective
import com.google.ortools.linearsolver.MPVariable
import org.apache.commons.math3.primes.Primes
import java.util.HashSet
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.roundToInt
import kotlin.system.measureTimeMillis

// MCMC where each Markov state is associated with the vertex that minimises
// an objective.
class MinimisationMCMC {
    val variables: Array<MPVariable>
    val solver = MPSolver("MinimisationMCMCSolver", MPSolver.OptimizationProblemType.GLOP_LINEAR_PROGRAMMING)
    val objectiveState: HashSet<Int>
    val currentSample: HashMap<Int,Double>

    val nConstraints: Int
        get() = solver.numConstraints()

    val nVariables: Int
        get() = variables.size

    val objective: ObjectiveWrapper

    class ObjectiveWrapper(val mpObjective: MPObjective, val variables: Array<MPVariable>) {
        operator fun get(i: Int) = mpObjective.getCoefficient(variables[i])
        operator fun set(i: Int, v: Double) { mpObjective.setCoefficient(variables[i], v) }
        fun clear() {
            mpObjective.clear()
            mpObjective.setMinimization()
        }
    }

//    constructor(constraints: List<Constraint<Number>>, initialSolution: Map<Int,Number>) {
//        variables = solver.makeNumVarArray(constraints.numVars()+constraints.numSlacks(), 0.0, Double.POSITIVE_INFINITY)
//        solver.setConstraintsWithExplicitSlacks(variables, constraints)
//        objectiveState = randomObjectiveState(nConstraints - initialSolution.size)
//        currentSample = HashMap(initialSolution.size)
//        initialSolution.forEach { currentSample[it.key] = it.value.toDouble() }
//    }

    constructor(constraints: List<Constraint<Number>>) {
        variables = solver.makeNumVarArray(constraints.numVars()+constraints.numSlacks(), 0.0, Double.POSITIVE_INFINITY)
        objective = ObjectiveWrapper(solver.objective(), variables)
        solver.setConstraintsWithExplicitSlacks(variables, constraints)
        currentSample = HashMap()
        println("Finding initial solution...")
        doMinimisation()
        objectiveState = HashSet()
    }


    fun nextSample(): Map<Int,Double> {
        val outgoingVariable = currentSample.keys.random()
        randomiseObjectiveState()
        setPrimeObjective(outgoingVariable)
        doMinimisation()

        return currentSample
    }


    fun setPrimeObjective(outgoingVariable: Int) {
        objective.clear()
        // first set chosen basis
        println("Calculating sum of primes")
        val sumOfLogPrimes = objectiveState.indices.sumByDouble { nthLogPrime(it+1) }
        var n = 0
//        val sumOfLogPrimesSquared = sumOfLogPrimes * sumOfLogPrimes
//        for(basisVar in objectiveState) {
//            val logPrime = nthLogPrime(++n)
//            objective[basisVar] = 0.0
//        }

        // now set fallbacks
        println("setting fallbacks ObjectiveState size is ${objectiveState.size}")
        n = 0
        var nextFreeBasis: Int
        var logPrime: Double
        for(basisVar in objectiveState) {
            logPrime = nthLogPrime(++n)
            nextFreeBasis = (basisVar + 1).rem(nVariables)
//            while(objectiveState.contains(nextFreeBasis) || objective[nextFreeBasis] != 0.0) {
            while(objectiveState.contains(nextFreeBasis) || objective[nextFreeBasis] != 0.0) {
                nextFreeBasis = (nextFreeBasis + 1).rem(nVariables)
            }
            objective[nextFreeBasis] = logPrime / sumOfLogPrimes
        }

        // now set excluded bases
//        val bigNumber = fallbackMultiplier * sumOfLogPrimes * 1.001
        println("setting unit bases")
        variables.indices.asSequence()
            .filter { !currentSample.containsKey(it) && !objectiveState.contains(it) && objective[it] == 0.0 }
            .forEach {
                objective[it] = 1.0
            }

        // now set outgoing variable
        objective[outgoingVariable] = 1.0

        println("Degeneracy is ${nConstraints - currentSample.size}")
        println("Sum of log primes is $sumOfLogPrimes")
        println("Smallest coefficient is ${ln(2.0)/sumOfLogPrimes}")
//        println("Big number is $bigNumber")
//        println("Set objective to ")
//        for(i in 0 until nVariables) {
//            print("${objective[i]} ")
//        }
        println()
    }


    fun setBinaryObjective(outgoingVariable: Int) {
        objective.clear()
        variables.indices.asSequence()
            .filter { !currentSample.contains(it) && !objectiveState.contains(it) && !objectiveState.contains(it-1) }
            .forEach {
                objective[it] = 1.0
            }

        // now set outgoing variable
        objective[outgoingVariable] = 1.0
    }


    fun doMinimisation() {
        var solveState: MPSolver.ResultStatus
        val solveTime = measureTimeMillis {
            solveState = solver.solve()
        }
        println("Solved in ${solveTime}ms")
        if (solveState == MPSolver.ResultStatus.OPTIMAL) {
            currentSample.clear()
            for(i in variables.indices) {
                val solution = variables[i].solutionValue()
                if(solution != 0.0) currentSample[i] = solution
            }
        } else
            throw(RuntimeException(solveState.name))
    }


    fun randomiseObjectiveState() {
        val degeneracy = nConstraints - currentSample.size
//        val unchosenVariables = variables.indices.toHashSet()
//        unchosenVariables.removeAll(currentSample.keys)
        println("Choosing random degeneracy state")
        objectiveState.clear()
        for(v in 1..(degeneracy+1)) {
            var nextChoice: Int
            do {
                nextChoice = variables.indices.random()
            } while(currentSample.containsKey(nextChoice) || objectiveState.contains(nextChoice))
            objectiveState.add(nextChoice)
//            unchosenVariables.remove(nextChoice)
        }
        println("Done")
    }

    companion object {

        init {
            System.loadLibrary("jniortools")
        }


        val logPrimes = arrayListOf(0.0, ln(2.0), ln(3.0))

        fun nthLogPrime(n: Int): Double {
            while(n >= logPrimes.size) {
                logPrimes.add(
                    ln(Primes
                        .nextPrime(exp(logPrimes.last()).roundToInt() + 1)
                        .toDouble())
                )
            }
            return logPrimes[n]
        }
    }


}