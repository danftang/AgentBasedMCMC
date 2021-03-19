package lib.sparseMatrix

import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector
import lib.sparseVector.asMutableVector
import java.lang.StringBuilder

class MutableSparseVectorGridMatrix<T>: MutableGridMatrix<T>, MutableEntryMatrix<T> {

    val _rows: ArrayList<MutableSparseVector<T>>
    val _columns: ArrayList<MutableSparseVector<T>>

    override val rows: List<SparseVector<T>>
        get() = _rows
    override val columns: List<SparseVector<T>>
        get() = _columns

    override val nonZeroEntries: Iterable<EntryMatrix.Entry<T>>
        get() = columns.asSequence()
            .withIndex()
            .flatMap { (j, col) ->
                col.nonZeroEntries.asSequence().map { (i, value) -> EntryMatrix.SimpleEntry(i,j,value) }
            }
            .asIterable()


    constructor(nRows: Int, nCols: Int, mutableSparseVectorFactory: () -> MutableSparseVector<T>) {
        this._rows = ArrayList(nRows)
        this._columns = ArrayList(nCols)
        repeat(nRows) { _rows.add(mutableSparseVectorFactory()) }
        repeat(nCols) { _columns.add(mutableSparseVectorFactory()) }
    }

    override fun set(row: Int, col: Int, value: T) {
        _rows[row][col] = value
        _columns[col][row] = value
    }

    override fun remove(row: Int, col: Int) {
        _rows[row].remove(col)
        _columns[col].remove(row)
    }

    override fun clear() {
        for(rowVector in _rows) rowVector.clear()
        for(colVector in _columns) colVector.clear()
    }

    override fun get(row: Int, col: Int): T {
        return columns[col][row]
    }

    override fun swapRows(row1: Int, row2: Int) {
        for(entry in rows[row1].nonZeroEntries) {
            _columns[entry.key][row2] = entry.value
            _columns[entry.key].remove(row1)
        }
        for(entry in rows[row2].nonZeroEntries) {
            _columns[entry.key][row1] = entry.value
            if(!_rows[row1].nonZeroEntries.containsKey(entry.key)) _columns[entry.key].remove(row2)
        }
        val tmpRow = _rows[row1]
        _rows[row1] = _rows[row2]
        _rows[row2] = tmpRow
    }

    override fun toString(): String {
        val out = StringBuilder()
        val colWidths = Array(nCols) { -1 }
        val elementStrings = Array(nRows) { row ->
            Array(nCols) { col ->
                val Mij = get(row,col)
                if(Mij != null) {
                    val elementStr = get(row, col).toString()
                    if (elementStr.length > colWidths[col])
                        colWidths[col] = elementStr.length
                    elementStr
                } else ""
            }
        }
        for(row in 0 until nRows) {
            for(col in 0 until nCols) {
                out.append(elementStrings[row][col].padEnd(colWidths[col],' '))
                out.append(" |")
            }
            out.append('\n')
        }
        return out.toString()
    }


}