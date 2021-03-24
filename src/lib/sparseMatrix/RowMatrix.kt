package lib.sparseMatrix

import lib.sparseVector.SparseVector

interface RowMatrix<T>: Matrix<T> {
    val rows: List<SparseVector<T>>


    // Default implementations
    /////////////////////////////////////////

    override val nRows: Int
        get() = rows.size
//
//    override operator fun get(row: Int, col: Int): T {
//        return rows[row][col]
//    }

}


operator fun <T> RowMatrix<T>.times(X: SparseVector<T>): SparseVector<T> {
    val result = X.new()
    for(i in rows.indices) {
        result[i] = rows[i].dotProduct(X)
    }
    return result
}