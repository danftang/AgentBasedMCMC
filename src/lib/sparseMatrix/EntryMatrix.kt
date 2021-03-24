package lib.sparseMatrix

import lib.abstractAlgebra.FieldOperators
import lib.sparseVector.MutableMapVector
import lib.sparseVector.SparseVector
import kotlin.math.min
import kotlin.math.sign

interface EntryMatrix<T>: Matrix<T> {
    val nonZeroEntries: Iterable<Entry<T>>

    interface Entry<out T> {
        val row: Int
        val col: Int
        val value: T
    }

    class SimpleEntry<out T>(override val row: Int, override val col: Int, override val value: T) : Entry<T>
}


//fun <T : Any> EntryMatrix<T>.diagonal(): SparseVector<T> {
//    val diag = MutableMapVector(operators)
//    for(i in 0 until min(nRows,nCols)) {
//        diag[i] = this[i,i]
//    }
//    return diag
//}


fun<T: Number> EntryMatrix<T>.toSparsityString(): String {
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


//inline fun<T: Any, R: Any> EntryMatrix<T>.mapNonZeroEntriesTo(destination: EntryMatrix<R>, transform: (T)->R): EntryMatrix<R> {
//    for(entry in nonZeroEntries) destination[entry.row, entry.col] = transform(entry.value)
//    return destination
//}
//
//inline fun<T: Any> EntryMatrix<T>.copyTo(destination: EntryMatrix<T>) {
//    for(entry in nonZeroEntries) destination[entry.row, entry.col] = entry.value
//}

