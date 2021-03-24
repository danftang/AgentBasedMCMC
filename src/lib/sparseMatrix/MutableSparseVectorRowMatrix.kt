package lib.sparseMatrix

import lib.sparseVector.MutableSparseVector

//class MutableSparseVectorRowMatrix<T> : MutableRowMatrix<T> {
//    override val nCols: Int
//
//    override val rows: ArrayList<MutableSparseVector<T>>
//
//    override val nRows: Int
//        get() = rows.size
//
//    constructor(nRows: Int, nCols: Int, mutableSparseVectorFactory: () -> MutableSparseVector<T>) {
//        this.nCols = nCols
//        this.rows = ArrayList<MutableSparseVector<T>>(nRows)
//        repeat(nRows) { rows.add(mutableSparseVectorFactory()) }
//    }
//
//    override fun set(row: Int, col: Int, value: T) {
//        rows[row][col] = value
//    }
//
//    override fun setToZero() {
//        for(rowVector in rows) rowVector.setToZero()
//    }
//
//    override fun get(row: Int, col: Int): T {
//        return rows[row][col]
//    }
//
//}