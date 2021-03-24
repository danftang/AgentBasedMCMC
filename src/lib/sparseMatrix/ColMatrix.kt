package lib.sparseMatrix

import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector

interface ColMatrix<T>: Matrix<T> {
    val columns: List<SparseVector<T>>



    // Default implementations
    /////////////////////////////////////////

    override val nCols: Int
        get() = columns.size

//    override operator fun get(row: Int, col: Int): T {
//        return columns[col][row]
//    }

}

fun <T> ColMatrix<T>.sparsity(): Double {
    val nonZeroEntriesSize = columns.sumBy { it.nonZeroEntries.size }
    return nonZeroEntriesSize.toDouble()/(nRows*nCols)
}
//operator fun <T> ColMatrix<T>.times(X: SparseVector<T>): SparseVector<T> {
//    val result =  X.new()
//    for(j in 0 until nCols) {
//        result.weightedPlusAssign(columns[j], X[j])
//    }
//    return result
//}