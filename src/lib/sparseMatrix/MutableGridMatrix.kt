package lib.sparseMatrix

import lib.sparseVector.SparseVector

interface MutableGridMatrix<T>: MutableRowMatrix<T>, MutableColMatrix<T> {
    fun swapRows(row1: Int, row2: Int)
}