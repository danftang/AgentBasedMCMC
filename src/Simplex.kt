import lib.sparseMatrix.ColMatrix
import lib.sparseMatrix.MutableColMatrix
import lib.sparseMatrix.MutableRowMatrix
import lib.sparseMatrix.RowMatrix
import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector

interface Simplex<T> {
    val M: Tableaux<T>
    val basicColsByRow: IntArray
    val B: MutableSparseVector<T>

    val objectiveRow: Int
        get() = M.nRows - 1
    val bColumn: Int
        get() = M.nCols - 1

    fun X(includeSlacks: Boolean): SparseVector<T>
    fun pivot(pivot: PivotPoint)
    fun pivotableRows(j: Int, allowPivotsOnNegativeElements: Boolean): List<Int>

}