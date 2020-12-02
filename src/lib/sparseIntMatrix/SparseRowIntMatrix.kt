package lib.sparseIntMatrix

interface SparseRowIntMatrix: SparseIntMatrix {
    val rows:       List<SparseIntVector>

    override val entries: Sequence<SparseIntMatrix.Entry>
        get() = rows
            .asSequence()
            .mapIndexed { i, col ->
                col.asSequence().map { SparseIntMatrix.Entry(i, it.key, it.value) }
            }
            .flatten()


    fun clearRow(i: Int)
    fun swapRows(i1: Int, i2: Int)

    fun removeRows(rowsToRemove: Iterable<Int>)

    fun replaceNonZeroElementsInRow(col: Int, map: (row: Int, value: Int) -> Int)

    // row[i1] = row[i1] + weight*row[i2]
    fun weightedRowPlusAssign(i1: Int, i2: Int, weight: Int) {
        for(entry in rows[i2]) {
            plusAssign(i1, entry.key, weight*entry.value)
        }
    }

    override operator fun times(X: SparseIntVector): HashIntVector {
        val result = HashIntVector()
        for(i in rows.indices) {
               result[i] = rows[i].dotProd(X)
        }
        return result
    }

    override operator fun times(X: IntArrayVector): IntArrayVector {
        return IntArrayVector(this.nRows) { i ->
            rows[i].dotProd(X)
        }
    }

    // Minimise CX
    // Subject to:
    // MX >= B
    // and
    // X >= 0 and X is integer
    // where M is this matrix
    // returns X
    //
//    fun IPsolve(B: SparseIntVector, C: List<Double> = DoubleArray(nCols) {0.0}.asList(), constraintType: String = ">="): IntArray {
//        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
//        val X = solver.makeIntVarArray(nCols, 0.0, Double.POSITIVE_INFINITY)
//        val constraints = Array<MPConstraint>(nRows) { solver.makeConstraint() }
//        for(col in 0 until nCols) {
//            for(coefficient in columns[col]) {
//                constraints[coefficient.key].setCoefficient(X[col], coefficient.value.toDouble())
//            }
//        }
//        for(i in 0 until nRows) {
//            when(constraintType) {
//                ">=" -> constraints[i].setBounds(B[i].toDouble(), Double.POSITIVE_INFINITY)
//                "=","==" -> constraints[i].setBounds(B[i].toDouble(), B[i].toDouble())
//                "<=" -> constraints[i].setBounds(Double.NEGATIVE_INFINITY, B[i].toDouble())
//                else -> throw(IllegalArgumentException("Unknown constraint type"))
//            }
//        }
//        val objective = solver.objective()
//        for(i in C.indices) {
//            objective.setCoefficient(X[i], C[i])
//        }
//        solver.objective().setMinimization()
//        val solveState: MPSolver.ResultStatus
//        val solveTime = measureTimeMillis {
//            solveState = solver.solve()
//        }
//        println("Solved in ${solveTime}ms")
//        return if (solveState == MPSolver.ResultStatus.OPTIMAL)
//            IntArray(X.size) { i -> X[i].solutionValue().toInt() }
//        else
//            throw(RuntimeException(
//                when (solveState) {
//                    MPSolver.ResultStatus.INFEASIBLE -> "Infeasible"
//                    MPSolver.ResultStatus.UNBOUNDED -> "Unbounded"
//                    else -> "Solve Error"
//                }
//            ))
//    }
//
//


}