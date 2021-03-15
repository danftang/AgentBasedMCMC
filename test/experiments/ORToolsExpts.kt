package experiments

import MutableConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPSolverParameters
import org.junit.Test
import kotlin.random.Random
import kotlin.system.measureTimeMillis

class ORToolsExpts {
    @Test
    fun solveTest() {
        val result = ORTools.GlopSolve(
            listOf(
                MutableConstraint(mutableMapOf(0 to 1),"<=", 1),
                MutableConstraint(mutableMapOf(0 to 1, 1 to 1),"<=", 2)
            ),
            mapOf(1 to -1)
        )
        println(result.asList())

        val intResult = ORTools.IntegerSolve(
            listOf(
                MutableConstraint(mutableMapOf(0 to 1),"<=", 1),
                MutableConstraint(mutableMapOf(0 to 1, 1 to 1),"<=", 2)
            ),
            mapOf(1 to -1)
        )
        println(intResult.asList())

    }

    @Test
    fun minTest() {
        System.loadLibrary("jniortools")
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.SCIP_MIXED_INTEGER_PROGRAMMING)
        val params = MPSolverParameters()
        params.setIntegerParam(MPSolverParameters.IntegerParam.INCREMENTALITY,
            MPSolverParameters.IncrementalityValues.INCREMENTALITY_ON.swigValue())
        println(MPSolverParameters.IncrementalityValues.INCREMENTALITY_ON.swigValue())
        println(MPSolverParameters.IncrementalityValues.INCREMENTALITY_OFF.swigValue())
        println(params.getIntegerParam(MPSolverParameters.IntegerParam.INCREMENTALITY))
//        solver.solve(params)
    }

    // Can we get incremental solving going?
    @Test
    fun incrementalSolveTest() {
        System.loadLibrary("jniortools")
        val nVars = 6000
        val nConstraints = 2000
        val solver = MPSolver("incremental?",MPSolver.OptimizationProblemType.SCIP_MIXED_INTEGER_PROGRAMMING)

        val vars = solver.makeIntVarArray(nVars, 0.0, 10.0)
        for(i in 1..nConstraints) {
            val mpConstraint = solver.makeConstraint()
            for(c in 0 until nVars) {
                mpConstraint.setCoefficient(vars[c], Random.nextInt(2).toDouble())
            }
        }
        val objective = solver.objective()
        for(c in 0 until nVars) {
            objective.setCoefficient(vars[c], Random.nextInt(2).toDouble())
        }
        objective.setMinimization()

        println("Starting Solve")
        // ######## FIRST SOLVE #########
        var solveState: MPSolver.ResultStatus
        val solveTime = measureTimeMillis {
            solveState = solver.solve()
        }
        println("${solveState.name} Solved in ${solveTime}ms")

//        val solution = DoubleArray(nVars) { i -> vars[i].solutionValue() }
//        solver.setHint(vars, solution)

        // ######### SECOND SOLVE ##########
        solver.objective().setCoefficient(vars[1], -1.0)
        val solveTime2 = measureTimeMillis {
            solveState = solver.solve()
        }
        println("${solveState.name} Solved in ${solveTime2}ms")


    }
}