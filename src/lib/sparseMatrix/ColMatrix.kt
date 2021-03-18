package lib.sparseMatrix

import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector

interface ColMatrix<T> {
    val columns: List<SparseVector<T>>
    val nRows: Int
    val nCols: Int
        get() = columns.size

}


operator fun <T> ColMatrix<T>.times(X: SparseVector<T>): SparseVector<T> {
    val result =  X.new()
    for(j in 0 until nCols) {
        result.weightedPlusAssign(columns[j], X[j])
    }
    return result
}