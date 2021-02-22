import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import kotlin.system.measureTimeMillis

object ORTools {
    init {
        System.loadLibrary("jniortools")
    }


    // Minimise CX
    // Subject to:
    // MX >= B
    // and
    // X >= 0 and X is integer
    // where M is this matrix
    // returns X
    //
    fun<T: Number> GlopSolve(constraints: List<Constraint<T>>, objective: Map<Int,T> = emptyMap()): DoubleArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.GLOP_LINEAR_PROGRAMMING)
        val nVariables = Integer.max(constraints.numVars(),
            objective.keys.max()?.let { it + 1 } ?: 0
        )
        val X = solver.makeNumVarArray(nVariables, 0.0, Double.POSITIVE_INFINITY)
        return ORSolve(solver, X, constraints, objective)
    }


    fun<T: Number> IntegerSolve(constraints: List<Constraint<T>>, objective: Map<Int,T> = emptyMap()): DoubleArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val nVariables = Integer.max(constraints.numVars(),
            objective.keys.max()?.let { it + 1 } ?: 0
        )
        val X = solver.makeIntVarArray(nVariables, 0.0, Double.POSITIVE_INFINITY)
        return ORSolve(solver, X, constraints, objective)
    }


    fun<T: Number> ORSolve(
        solver: MPSolver,
        X: Array<MPVariable>,
        constraints: List<Constraint<T>>,
        objective: Map<Int,T> = emptyMap()): DoubleArray {
        for(constraint in constraints) {
            val mpConstraint = solver.makeConstraint()
            for((varId, coeff) in constraint.coefficients) {
                mpConstraint.setCoefficient(X[varId], coeff.toDouble())
            }
            val b = constraint.constant.toDouble()
            when (constraint.relation) {
                ">=" -> mpConstraint.setBounds(b, Double.POSITIVE_INFINITY)
                "=", "==" -> mpConstraint.setBounds(b, b)
                "<=" -> mpConstraint.setBounds(Double.NEGATIVE_INFINITY, b)
                else -> throw(IllegalArgumentException("Unknown constraint type"))
            }
        }

        val mpObjective = solver.objective()
        for (entry in objective) {
            mpObjective.setCoefficient(X[entry.key], entry.value.toDouble())
        }
        mpObjective.setMinimization()
        var solveState: MPSolver.ResultStatus
        val solveTime = measureTimeMillis {
            solveState = solver.solve()
        }
        println("Solved in ${solveTime}ms")
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

}