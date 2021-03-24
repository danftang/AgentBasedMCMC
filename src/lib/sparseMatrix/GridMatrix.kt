package lib.sparseMatrix

import lib.abstractAlgebra.FieldOperators
import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector
import lib.sparseVector.asMutableVector
import java.lang.StringBuilder

class GridMatrix<T> private constructor(field: FieldOperators<T>):
    MutableGridMatrix<T>, MutableEntryMatrix<T>, FieldOperators<T> by field {

    private val _rows = ArrayList<MutableSparseVector<T>>()
    private val _columns = ArrayList<MutableSparseVector<T>>()

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


    constructor(nRows: Int, nCols: Int, mutableSparseVectorFactory: () -> MutableSparseVector<T>):
            this(mutableSparseVectorFactory().operators) {
        _rows.ensureCapacity(nRows)
        _columns.ensureCapacity(nCols)
        for(i in 0 until nRows) {
            _rows.add(mutableSparseVectorFactory())
        }
        for(j in 0 until nCols) {
            _columns.add(mutableSparseVectorFactory())
        }
    }


    constructor(nRows: Int, nCols: Int, operators: FieldOperators<T>): this(nRows, nCols, {
        HashMap<Int,T>().asMutableVector(operators)
    })


    override fun set(row: Int, col: Int, value: T) {
        _rows[row][col] = value
        _columns[col][row] = value
    }

    override fun setToZero() {
        for(rowVector in _rows) rowVector.setToZero()
        for(colVector in _columns) colVector.setToZero()
    }

    override fun get(row: Int, col: Int): T {
        return columns[col][row]
    }

    override fun setRowToZero(row: Int) {
        for(entry in _rows[row].nonZeroEntries) {
            _columns[entry.key][row] = zero
        }
        _rows[row].setToZero()
    }

    override fun mapAssignRow(row: Int, transform: (T) -> T) {
        _rows[row].mapAssignWithIndex { col, Mij ->
            val newMij = transform(Mij)
            _columns[col][row] = newMij
            newMij
        }
    }

    override fun mapAssignRowWithIndex(row: Int, transform: (Int, T) -> T) {
        _rows[row].mapAssignWithIndex { col, Mij ->
            val newMij = transform(col, Mij)
            _columns[col][row] = newMij
            newMij
        }
    }

    override fun setColToZero(col: Int) {
        for(entry in _columns[col].nonZeroEntries) {
            _rows[entry.key][col] = zero
        }
        _columns[col].setToZero()
    }

    override fun mapAssignCol(col: Int, transform: (T) -> T) {
        _columns[col].mapAssignWithIndex { row, Mij ->
            val newMij = transform(Mij)
            _rows[row][col] = newMij
            newMij
        }
    }

    override fun mapAssignColWithIndex(col: Int, transform: (Int, T) -> T) {
        _columns[col].mapAssignWithIndex { row, Mij ->
            val newMij = transform(row,Mij)
            _rows[row][col] = newMij
            newMij
        }
    }


    override fun swapRows(row1: Int, row2: Int) {
        for(entry in rows[row1].nonZeroEntries) {
            _columns[entry.key][row2] = entry.value
            _columns[entry.key][row1] = zero
        }
        for(entry in rows[row2].nonZeroEntries) {
            _columns[entry.key][row1] = entry.value
            if(!_rows[row1].nonZeroEntries.containsKey(entry.key)) _columns[entry.key][row2] = zero
        }
        val tmpRow = _rows[row1]
        _rows[row1] = _rows[row2]
        _rows[row2] = tmpRow
    }

    fun getMutableRowView(row: Int): MutableSparseVector<T> = RowView(row)

    fun getMutableColView(col: Int): MutableSparseVector<T> = ColView(col)

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

    inner class RowView(val row: Int): MutableSparseVector<T>, FieldOperators<T> by operators {
        override val nonZeroEntries: Map<Int, T>
            get() = _rows[row].nonZeroEntries

        override fun set(index: Int, value: T) { this@GridMatrix[row,index] = value }
        override fun setToZero() { setRowToZero(row) }
        override fun new(): MutableSparseVector<T> { return _rows[row].new() }
        override fun mapAssign(elementTransform: (T) -> T) { mapAssignRow(row,elementTransform) }
        override fun mapAssignWithIndex(transform: (Int, T) -> T) { mapAssignRowWithIndex(row, transform) }
    }

    inner class ColView(val col: Int): MutableSparseVector<T>, FieldOperators<T> by operators {
        override val nonZeroEntries: Map<Int, T>
            get() = _columns[col].nonZeroEntries

        override fun set(index: Int, value: T) { this@GridMatrix[index,col] = value }
        override fun setToZero() { setColToZero(col) }
        override fun new(): MutableSparseVector<T> { return _columns[col].new() }
        override fun mapAssign(elementTransform: (T) -> T) { mapAssignCol(col, elementTransform) }
        override fun mapAssignWithIndex(transform: (Int, T) -> T) { mapAssignColWithIndex(col, transform) }
    }
}