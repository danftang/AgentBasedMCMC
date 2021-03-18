import lib.sparseMatrix.ColMatrix
import lib.sparseMatrix.MutableColMatrix
import lib.sparseMatrix.MutableRowMatrix
import lib.sparseMatrix.RowMatrix
import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector

interface Simplex<T,MATRIX>: RowSimplex<T> where MATRIX: MutableRowMatrix<T>, MATRIX: ColMatrix<T> {
    override val M: MATRIX
}