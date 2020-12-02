package lib.sparseMatrix

import lib.vector.SparseVector

interface SparseRowColMatrix<T: Any>: SparseRowMatrix<T>, SparseColMatrix<T> {
    override val entries: Sequence<SparseMatrix.Entry<T>>
        get() = (this as SparseColMatrix<T>).entries

    override fun get(row: Int, col: Int)        = super<SparseColMatrix>.get(row,col)
    override fun times(X: SparseVector<T>)      = super<SparseColMatrix>.times(X)
    override fun set(row: Int, col: Int, value: T) = super<SparseColMatrix>.set(row,col,value)
    override fun mapAssign(row: Int, col: Int, remappingFunction: (T) -> T) = super<SparseColMatrix>.mapAssign(row,col,remappingFunction)
}