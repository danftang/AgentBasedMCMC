package lib.sparseMatrix

import lib.abstractAlgebra.FieldOperators
import lib.collections.GridMap
import lib.sparseVector.MutableMapVector
import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector

class GridMapMatrix<T: Any>(override val operators: FieldOperators<T>, val gridMap: GridMap<T>): SparseRowMatrix<T>, SparseColMatrix<T>, FieldOperators<T> by operators {
    override val rows = VectorizedList(gridMap.rows, operators)
    override val columns = VectorizedList(gridMap.columns, operators)
    override val nonZeroEntries: Iterable<SparseMatrix.Entry<T>>
        get() = super<SparseColMatrix>.nonZeroEntries

    constructor(operators: FieldOperators<T>, nRows: Int, nCols: Int): this(operators, GridMap(nRows, nCols))

    override fun get(row: Int, col: Int)        = gridMap.get(row,col)?:zero
    override fun times(X: SparseVector<T>)      = super<SparseColMatrix>.times(X)
    override fun set(row: Int, col: Int, value: T) {
        if(value == zero) gridMap.remove(row,col) else gridMap.set(row,col,value)
    }
    override fun mapAssign(row: Int, col: Int, remappingFunction: (T) -> T) = super<SparseColMatrix>.mapAssign(row,col,remappingFunction)

    fun resize(nRows: Int, nCols: Int) { gridMap.resize(nRows,nCols) }

    override fun toString(): String {
        return gridMap.toString()
    }

    class VectorizedList<T>(val list: List<MutableMap<Int,T>>, val operators: FieldOperators<T>): AbstractList<MutableSparseVector<T>>() {
        override val size: Int
            get() = list.size
        override fun get(index: Int) = MutableMapVector(operators, list[index])
    }

}