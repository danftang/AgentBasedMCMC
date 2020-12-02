package lib.sparseMatrix

import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import lib.abstractAlgebra.FieldOperators
import lib.vector.MutableMapVector
import lib.vector.SparseVector
import java.lang.IllegalArgumentException
import kotlin.math.min
import kotlin.math.sign
import kotlin.system.measureTimeMillis

interface SparseMatrix<T: Any>: FieldOperators<T> {
    data class Entry<T>(val row: Int, val col: Int, val value: T)

    val entries: Sequence<Entry<T>>

    val nRows: Int
    val nCols: Int

    operator fun get(row: Int, col: Int): T
    operator fun times(X: SparseVector<T>): SparseVector<T>

    operator fun set(row: Int, col: Int, value: T)

    // remappingFunction takes the old value and returns the new
    fun mapAssign(row: Int, col: Int, remappingFunction: (T)->T)

//    fun newSparseVector(): MutableSparseVector<T>
//    fun MutableMap<Int,T>.asMutableSparseVector(): MutableSparseVector<T>

    fun diagonal(): SparseVector<T> {
        val diag = MutableMapVector(this)
        for(i in 0 until min(nRows,nCols)) {
            diag[i] = this[i,i]
        }
        return diag
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
fun<T: Number> SparseMatrix<T>.IPsolve(B: SparseVector<T>, C: List<Double> = DoubleArray(nCols) {0.0}.asList(), constraintType: String = ">="): DoubleArray {
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



fun<T: Number> SparseMatrix<T>.toSparsityString(): String {
    val out = StringBuilder()
    for(row in 0 until nRows) {
        for (col in 0 until nCols) {
            val v = this[row,col]
            out.append(when(v.toDouble().sign) {
                1.0 -> '+'
                -1.0 -> '-'
                else ->  '.'
            })
        }
        out.appendln()
    }
    return out.toString()
}

