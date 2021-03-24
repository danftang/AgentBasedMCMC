import lib.sparseMatrix.MutableColMatrix
import lib.sparseMatrix.MutableRowMatrix
import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector

interface RowSimplex<T> {
    val basicColsByRow: IntArray
    val M: MutableRowMatrix<T>
    val B: MutableSparseVector<T>

    val objectiveRow: Int
        get() = M.nRows - 1
    val bColumn: Int
        get() = M.nCols - 1

    fun X(includeSlacks: Boolean): SparseVector<T>
    fun pivot(pivot: PivotPoint)
    fun pivotableRows(j: Int, allowPivotsOnNegativeElements: Boolean): List<Int>

}

