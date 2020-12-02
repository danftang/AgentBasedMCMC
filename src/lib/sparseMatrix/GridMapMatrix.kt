package lib.sparseMatrix

import lib.abstractAlgebra.FieldOperators
import lib.collections.MutableViewGridMap
import lib.vector.MutableMapVector
import lib.vector.MutableSparseVector

class GridMapMatrix<T: Any>(val operators: FieldOperators<T>, val gridMap: MutableViewGridMap<T>): SparseRowColMatrix<T>, FieldOperators<T> by operators {
    override val rows = VectorizedList(gridMap.rows)
    override val columns = VectorizedList(gridMap.columns)

    constructor(operators: FieldOperators<T>, nRows: Int, nCols: Int): this(operators, MutableViewGridMap(nRows, nCols))

    fun resize(nRows: Int, nCols: Int) { gridMap.resize(nRows,nCols) }

    override fun toString(): String {
        return gridMap.toString()
    }

    inner class VectorizedList(val list: List<MutableMap<Int,T>>): AbstractList<MutableSparseVector<T>>() {
        override val size: Int
            get() = list.size
        override fun get(index: Int) = MutableMapVector(this@GridMapMatrix, list[index])
    }
}