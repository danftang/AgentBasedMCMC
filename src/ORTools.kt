import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import lib.sparseMatrix.GridMapMatrix
import lib.sparseMatrix.SparseMatrix
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


    fun<T: Number> GlopSolve(tableaux: SparseMatrix<T>): DoubleArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.GLOP_LINEAR_PROGRAMMING)
        val bColumn = tableaux.nCols - 1
        val objectiveRow = tableaux.nRows - 1
        val X = solver.makeNumVarArray(tableaux.nCols - 1, 0.0, Double.POSITIVE_INFINITY)
        val constraints = Array<MPConstraint>(tableaux.nRows-1) { solver.makeConstraint() }
        val objective = solver.objective()
        for(entry in tableaux.nonZeroEntries) {
            if(entry.row != objectiveRow) {
                if(entry.col != bColumn) {
                    constraints[entry.row].setCoefficient(X[entry.col], entry.value.toDouble())
                } else {
                    constraints[entry.row].setBounds(entry.value.toDouble(), entry.value.toDouble())
                }
            } else {
                objective.setCoefficient(X[entry.col], entry.value.toDouble())
            }
        }
        solver.objective().setMinimization()
        return doSolve(solver,X)
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
        return doSolve(solver, X)
    }


    fun doSolve(solver: MPSolver, X: Array<MPVariable>): DoubleArray {
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