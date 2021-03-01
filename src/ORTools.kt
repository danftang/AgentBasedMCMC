import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import lib.abstractAlgebra.DoubleOperators
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
        val X = solver.makeNumVarArray(tableaux.nCols - 1, 0.0, Double.POSITIVE_INFINITY)
        return ORSolve(solver, X, tableaux)
    }


    fun<T: Number> IntegerSolve(constraints: List<Constraint<T>>, objective: Map<Int,T> = emptyMap()): DoubleArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val nVariables = Integer.max(constraints.numVars(),
            objective.keys.max()?.let { it + 1 } ?: 0
        )
        val X = solver.makeIntVarArray(nVariables, 0.0, Double.POSITIVE_INFINITY)
        return ORSolve(solver, X, constraints, objective)
    }


    fun<T: Number> IntegerSolve(tableaux: SparseMatrix<T>): DoubleArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val X = solver.makeIntVarArray(tableaux.nCols-1, 0.0, Double.POSITIVE_INFINITY)
        return ORSolve(solver, X, tableaux)
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


    fun<T: Number> ORSolve(
        solver: MPSolver,
        X: Array<MPVariable>,
        tableaux: SparseMatrix<T>): DoubleArray {
        val bColumn = tableaux.nCols - 1
        val objectiveRow = tableaux.nRows - 1
        val constraints = Array<MPConstraint>(tableaux.nRows-1) { solver.makeConstraint(0.0,0.0) }
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


    fun MPSolver.toSimplex(): Simplex<Double> {
        val grid = GridMapMatrix<Double>(DoubleOperators,numConstraints()+1, numVariables()+1)

        // check format
        constraints().forEach { constraint ->
            if(constraint.lb() != constraint.ub()) throw(java.lang.RuntimeException("Can only do equality constraints at the moment"))
        }
        variables().forEach { variable ->
            if(variable.lb() != 0.0 || variable.ub() != Double.POSITIVE_INFINITY) throw(java.lang.RuntimeException("Can only do positive variables for now"))
        }

        constraints().forEach { constraint ->
            variables().forEach { variable ->
                val coeff = constraint.getCoefficient(variable)
                if(coeff != 0.0) grid[constraint.index(), variable.index()] = coeff
            }
            grid[constraint.index(), grid.nCols-1] = constraint.lb()
        }

        variables().forEach { variable ->
            val coeff = objective().getCoefficient(variable)
            grid[grid.nRows-1, variable.index()] = coeff
        }

        val simplex = Simplex(grid)

//    this.solve()
//    val initialSolution = variables()
//        .mapNotNull { if(it.solutionValue() != 0.0) it.index() else null }
//
//    simplex.pivotColumns(initialSolution)

        return simplex
    }

}