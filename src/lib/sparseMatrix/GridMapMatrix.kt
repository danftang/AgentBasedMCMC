package lib.sparseMatrix

import lib.abstractAlgebra.FieldOperators
import lib.collections.GridMap
import lib.sparseVector.MutableMapVector
import lib.sparseVector.SparseVector

class GridMapMatrix<T: Any>(override val operators: FieldOperators<T>, val gridMap: GridMap<T>): MutableGridMatrix<T>, EntryMatrix<T>, FieldOperators<T> by operators {

    override val rows = VectorizedList(gridMap.rows)

    override val columns = VectorizedList(gridMap.columns)

    override val nonZeroEntries: Iterable<EntryMatrix.Entry<T>>
        get() = gridMap.entries.asSequence().map { Entry(it) }.asIterable()

    override val nRows: Int
        get() = gridMap.nRows

    override val nCols: Int
        get() = gridMap.nCols

    constructor(operators: FieldOperators<T>, nRows: Int, nCols: Int): this(operators, GridMap(nRows, nCols))

    override operator fun get(row: Int, col: Int)        = gridMap[row,col]?:zero


    override operator fun set(row: Int, col: Int, value: T) {
        if(value.isZero()) gridMap.remove(row,col) else gridMap[row,col] = value
    }

    override fun remove(row: Int, col: Int) {
        gridMap.remove(row,col)
    }

    override fun clear() {
        gridMap.clear()
    }


    fun resize(nRows: Int, nCols: Int) { gridMap.resize(nRows,nCols) }

//    operator fun times(X: SparseVector<T>): SparseVector<T>      = (this as ColMatrix<T>) * X
//    override fun mapAssign(row: Int, col: Int, remappingFunction: (T) -> T) = super<MutableColMatrix>.mapAssign(row,col,remappingFunction)

    override fun toString(): String {
        return gridMap.toString()
    }

    inner class VectorizedList(val list: List<MutableMap<Int,T>>): AbstractList<SparseVector<T>>() {
        override val size: Int
            get() = list.size
        override fun get(index: Int) = MutableMapVector(operators, list[index])
    }

    inline class Entry<T>(val gridMapEntry: GridMap.GridMapEntry<T>): EntryMatrix.Entry<T> {
        override val row: Int
            get() = gridMapEntry.row
        override val col: Int
            get() = gridMapEntry.col
        override val value: T
            get() = gridMapEntry.value
    }

    override fun swapRows(row1: Int, row2: Int) {
        gridMap.swapRows(row1, row2)
    }


}