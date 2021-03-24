package lib.sparseMatrix

interface Matrix<T> {
    val nRows: Int
    val nCols: Int

    operator fun get(row: Int, col: Int): T
}