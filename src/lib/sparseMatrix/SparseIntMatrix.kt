package lib.sparseMatrix

import kotlin.math.min
import kotlin.math.sign

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
    operator fun times(X: IntVector): IntVector

    fun diagonal(): HashIntVector {
        val diag = HashIntVector()
        for(i in 0 until min(nRows,nCols)) {
            diag[i] = this[i,i]
        }
        return diag
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