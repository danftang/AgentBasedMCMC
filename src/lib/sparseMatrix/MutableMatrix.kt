package lib.sparseMatrix

interface MutableMatrix<T>: Matrix<T> {
    operator fun set(row: Int, col: Int, value: T)
    fun remove(row: Int, col: Int)
    fun clear()

}