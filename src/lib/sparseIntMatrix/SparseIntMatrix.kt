package lib.sparseIntMatrix

import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import java.lang.IllegalArgumentException
import kotlin.math.min
import kotlin.math.sign
import kotlin.system.measureTimeMillis

interface SparseIntMatrix {
    data class Entry(val row: Int, val col: Int, val value: Int)

    val entries: Sequence<Entry>

    val nRows: Int
    val nCols: Int

    operator fun set(row: Int, col: Int, value: Int)
    operator fun get(row: Int, col: Int): Int
    fun plusAssign(row: Int, col: Int, addition: Int)
    fun resize(nRows: Int, nCols: Int)

    operator fun times(X: SparseIntVector): SparseIntVector
    operator fun times(X: IntArrayVector): IntArrayVector

    fun diagonal(): HashIntVector {
        val diag = HashIntVector()
        for(i in 0 until min(nRows,nCols)) {
            diag[i] = this[i,i]
        }
        return diag
    }

    // Minimise CX
    // Subject to:
    // MX >= B
    // and
    // X >= 0 and X is integer
    // where M is this matrix
    // returns X
    //
    fun IPsolve(B: SparseIntVector, C: List<Double> = DoubleArray(nCols) {0.0}.asList(), constraintType: String = ">="): IntArrayVector {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val X = solver.makeIntVarArray(nCols, 0.0, Double.POSITIVE_INFINITY)
        val constraints = Array<MPConstraint>(nRows) { solver.makeConstraint() }
        for(entry in entries) {
            constraints[entry.row].setCoefficient(X[entry.col], entry.value.toDouble())
        }
        for(i in 0 until nRows) {
            when(constraintType) {
                ">=" -> constraints[i].setBounds(B[i].toDouble(), Double.POSITIVE_INFINITY)
                "=","==" -> constraints[i].setBounds(B[i].toDouble(), B[i].toDouble())
                "<=" -> constraints[i].setBounds(Double.NEGATIVE_INFINITY, B[i].toDouble())
                else -> throw(IllegalArgumentException("Unknown constraint type"))
            }
        }
        val objective = solver.objective()
        for(i in C.indices) {
            objective.setCoefficient(X[i], C[i])
        }
        solver.objective().setMinimization()
        var solveState: MPSolver.ResultStatus? = null
        val solveTime = measureTimeMillis {
            solveState = solver.solve()
        }
        println("Solved in ${solveTime}ms")
        return if (solveState == MPSolver.ResultStatus.OPTIMAL)
            IntArrayVector(X.size) { i -> X[i].solutionValue().toInt() }
        else
            throw(RuntimeException(
                when (solveState) {
                    MPSolver.ResultStatus.INFEASIBLE -> "Infeasible"
                    MPSolver.ResultStatus.UNBOUNDED -> "Unbounded"
                    else -> "Solve Error"
                }
            ))
    }


    fun toSparsityString(): String {
        val out = StringBuilder()
        for(row in 0 until nRows) {
            for (col in 0 until nCols) {
                val v = this[row,col]
                out.append(when(v.sign) {
                    1 -> '+'
                    -1 -> '-'
                    else ->  '.'
                })
            }
            out.appendln()
        }
        return out.toString()
    }


}