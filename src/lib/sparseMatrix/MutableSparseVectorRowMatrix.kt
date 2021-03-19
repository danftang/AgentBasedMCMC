package lib.sparseMatrix

import lib.abstractAlgebra.FieldOperators
import lib.collections.IntToDoubleHashMap
import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector
import lib.sparseVector.asMutableVector
import org.apache.commons.math3.util.OpenIntToDoubleHashMap
import kotlin.math.absoluteValue

class MutableSparseVectorRowMatrix<T> : MutableRowMatrix<T> {
    override val nCols: Int

    override val rows: ArrayList<MutableSparseVector<T>>

    override val nRows: Int
        get() = rows.size

    constructor(nRows: Int, nCols: Int, mutableSparseVectorFactory: () -> MutableSparseVector<T>) {
        this.nCols = nCols
        this.rows = ArrayList<MutableSparseVector<T>>(nRows)
        repeat(nRows) { rows.add(mutableSparseVectorFactory()) }
    }

    override fun set(row: Int, col: Int, value: T) {
        rows[row][col] = value
    }

    override fun remove(row: Int, col: Int) {
        rows[row].remove(col)
    }

    override fun clear() {
        for(rowVector in rows) rowVector.clear()
    }

    override fun get(row: Int, col: Int): T {
        return rows[row][col]
    }

}