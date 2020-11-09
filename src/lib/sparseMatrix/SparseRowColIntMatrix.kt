package lib.sparseMatrix

import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import java.lang.IllegalArgumentException
import kotlin.math.min
import kotlin.math.sign
import kotlin.system.measureTimeMillis

//interface SparseRowColIntMatrix {
//
//    val columns:    List<SparseIntVector>
//    val rows:       List<SparseIntVector>
//
//    val nRows: Int
//    val nCols: Int
//
//    operator fun set(row: Int, col: Int, value: Int)
//    operator fun get(row: Int, col: Int): Int
//    fun plusAssign(row: Int, col: Int, addition: Int)
//    fun clearRow(i: Int)
//    fun clearColumn(j: Int)
//    fun swapRows(i1: Int, i2: Int)
//    fun resize(nRows: Int, nCols: Int)
//    fun copy(): SparseRowColIntMatrix
//
//    fun setColumn(col: Int, vector: Iterable<Map.Entry<Int,Int>>) {
//        clearColumn(col)
//        for(entry in vector) {
//            this[entry.key, col] = entry.value
//        }
//    }
//
//    fun weightedRowPlusAssign(i1: Int, i2: Int, weight: Int) {
//        for(entry in rows[i2]) {
//            plusAssign(i1, entry.key, weight*entry.value)
//        }
//    }
//
//    fun weightedColPlusAssign(i1: Int, i2: Int, weight: Int) {
//        for(entry in columns[i2]) {
//            plusAssign(entry.key, i1, weight*entry.value)
//        }
//    }
//
//
//    operator fun times(X: SparseIntVector): HashIntVector {
//        val result = HashIntVector()
//        for(entry in X) {
//            result.weightedPlusAssign(columns[entry.key], entry.value)
//        }
//        return result
//    }
//
//    operator fun times(X: IntArray): SparseIntVector {
//        val result = HashIntVector()
//        for(i in X.indices) {
//            if(X[i] != 0) result.weightedPlusAssign(columns[i], X[i])
//        }
//        return result
//    }
//
//    operator fun times(X: List<Int>): SparseIntVector {
//        val result = HashIntVector()
//        for(i in X.indices) {
//            if(X[i] != 0) result.weightedPlusAssign(columns[i], X[i])
//        }
//        return result
//    }
//
//
//    fun diagonal(): HashIntVector {
//        val diag = HashIntVector()
//        for(i in 0 until min(nRows,nCols)) {
//            diag[i] = this[i,i]
//        }
//        return diag
//    }
//
//    fun toSparsityString(): String {
//        val out = StringBuilder()
//        for(row in 0 until nRows) {
//            for (col in columns) {
//                val v = col[row]
//                out.append(when(v.sign) {
//                    1 -> '+'
//                    -1 -> '-'
//                    else ->  '.'
//                })
//            }
//            out.appendln()
//        }
//        return out.toString()
//    }
//
//    // Minimise CX
//    // Subject to:
//    // MX >= B
//    // and
//    // X >= 0 and X is integer
//    // where M is this matrix
//    // returns X
//    //
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
//}