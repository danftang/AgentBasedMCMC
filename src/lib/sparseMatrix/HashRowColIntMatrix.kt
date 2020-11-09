package lib.sparseMatrix

import java.util.AbstractMap

class HashRowColIntMatrix: SparseRowIntMatrix, SparseColIntMatrix {
    override val columns: List<SparseIntVector>
        get() = _columns
    override val rows: List<SparseIntVector>
        get() = _rows

    override val entries: Sequence<SparseIntMatrix.Entry>
        get() = super<SparseColIntMatrix>.entries

    override val nRows: Int
        get() = rows.size
    override val nCols: Int
        get() = columns.size

    private val _columns: ArrayList<Vector>
    private val _rows: ArrayList<Vector>


    constructor(nRows: Int, nCols: Int) {
        _columns = ArrayList(nCols)
        _rows = ArrayList(nRows)
        for(j in 1..nCols) _columns.add(Vector())
        for(i in 1..nRows) _rows.add(Vector())
    }


    constructor(copy: SparseIntMatrix): this(copy.nRows, copy.nCols) {
        for(entry in copy.entries) {
            this[entry.row, entry.col] = entry.value
        }
    }


    override operator fun set(row: Int, col: Int, value: Int) {
        if(value == 0) {
            remove(row,col)
            return
        }
        _rows[row].data.getOrPut(col) {
            val newEntry = EntryValue(value)
            _columns[col].data[row] = newEntry
            newEntry
        }.value = value
    }

    override operator fun get(row: Int, col: Int): Int = rows[row][col]

    private fun entry(row: Int, col: Int): EntryValue? = _rows[row].data[col]

    override fun plusAssign(row: Int, col: Int, addition: Int) {
        if(addition == 0) return
        entry(row,col)
            ?.run {
                value += addition
                if(value == 0) remove(row,col)
            }
            ?:run {
                val newEntry = EntryValue(addition)
                _rows[row][col] = newEntry
                _columns[col][row] = newEntry
            }
    }


    override fun clearRow(i: Int) {
        for(entry in rows[i]) _columns[entry.key].remove(i)
        _rows[i].data.clear()
    }

    override fun clearColumn(j: Int) {
        for(entry in columns[j]) _rows[entry.key].remove(j)
        _columns[j].data.clear()
    }

    override fun times(X: SparseIntVector) = super<SparseColIntMatrix>.times(X)

    override fun times(X: IntVector) = super<SparseColIntMatrix>.times(X)

    operator fun times(M: SparseIntMatrix): HashRowColIntMatrix {
        assert(nCols == M.nRows)
        val product = HashRowColIntMatrix(nRows, M.nCols)
        for(otherEntry in M.entries) {
            for (entry in columns[otherEntry.row]) {
                product[entry.key, otherEntry.col] += entry.value * otherEntry.value
            }
        }
        return product
    }

    override fun removeColumns(colsToRemove: Iterable<Int>) {
        _columns.removeAll(colsToRemove.map { _columns[it] } )
        _rows.forEach { it.data.clear() }
        for(thisCol in 0 until nCols) {
            for(entry in _columns[thisCol].data) {
                _rows[entry.key][thisCol] = entry.value
            }
        }
    }

    override fun swapRows(i1: Int, i2: Int) {
        val tmp = _rows[i1]
        _rows[i1] = _rows[i2]
        _rows[i2] = tmp
        for(entry in _rows[i1].data) _columns[entry.key].remove(i2)
        for(entry in _rows[i2].data) _columns[entry.key].remove(i1)
        for(entry in _rows[i1].data) _columns[entry.key][i1] = entry.value
        for(entry in _rows[i2].data) _columns[entry.key][i2] = entry.value
    }

    override fun swapCols(i1: Int, i2: Int) {
        val tmp = _columns[i1]
        _columns[i1] = _columns[i2]
        _columns[i2] = tmp
        for(entry in _columns[i1].data) _rows[entry.key].remove(i2)
        for(entry in _columns[i2].data) _rows[entry.key].remove(i1)
        for(entry in _columns[i1].data) _rows[entry.key][i1] = entry.value
        for(entry in _columns[i2].data) _rows[entry.key][i2] = entry.value
    }


    override fun removeRows(rowsToRemove: Iterable<Int>) {
        _rows.removeAll(rowsToRemove.map { _columns[it] } )
        _columns.forEach { it.data.clear() }
        for(thisRow in 0 until nRows) {
            for(entry in _rows[thisRow].data) {
                _columns[entry.key][thisRow] = entry.value
            }
        }
    }

    // resize this matrix, adding zeroes when expanding
    // and deleting the higher index rows/columns when
    // shrinking, while leaving other entries unchanged
    override fun resize(nRows: Int, nCols: Int) {
        if(nRows > rows.size) {
            repeat(nRows - rows.size) { _rows.add(Vector()) }
        } else if(nRows < rows.size) {
            repeat(rows.size - nRows) {
                clearRow(_rows.lastIndex)
                _rows.removeAt(_rows.lastIndex)
            }
        }
        if(nCols > columns.size) {
            repeat(nCols - columns.size) { _columns.add(Vector()) }
        } else if(nCols < columns.size) {
            repeat(columns.size - nCols) {
                clearColumn(_columns.lastIndex)
                _columns.removeAt(_columns.lastIndex)
            }
        }
    }


    override fun replaceNonZeroElementsInCol(col: Int, map: (row: Int, value: Int) -> Int) {
        _columns[col].data.forEach { entry ->
            entry.value.value = map(entry.key, entry.value.value)
        }
    }


    override fun replaceNonZeroElementsInRow(row: Int, map: (row: Int, value: Int) -> Int) {
        _rows[row].data.forEach { entry ->
            entry.value.value = map(entry.key, entry.value.value)
        }
    }


    fun timesAssignCol(col: Int, multiplier: Int) {
        _columns[col].data.forEach { entry ->
            entry.value.value = entry.value.value * multiplier
        }
    }

    fun timesAssignRow(row: Int, multiplier: Int) {
        _rows[row].data.forEach { entry ->
            entry.value.value = entry.value.value * multiplier
        }
    }


    private fun remove(row: Int, col: Int) {
        _rows[row].remove(col)
        _columns[col].remove(row)
    }


    fun copy() = HashRowColIntMatrix(this)

    override fun toString(): String {
        val out = StringBuilder()
        for(row in rows.indices) {
            for (col in columns.indices) {
                val v = this[row,col]
                if(v in 0..9) out.append(' ')
                out.append(v)
                out.append(' ')
            }
            out.append('\n')
        }
        return out.toString()
    }



    data class EntryValue(var value: Int)

    class Vector(val data: HashMap<Int, EntryValue> = HashMap(4)) : SparseIntVector {
        override val sparseSize: Int
            get() = data.size

        override operator fun get(index: Int): Int {
            return data[index]?.value?:0
        }

        override fun iterator(): Iterator<Map.Entry<Int, Int>> =
            data.entries.asSequence().map { AbstractMap.SimpleEntry(it.key, it.value.value) }.iterator()


        operator fun set(index: Int, entry: EntryValue) {
            data[index] = entry
        }

        fun remove(index: Int) {
            data.remove(index)
        }

        override fun toString(): String {
            return data.mapValues { it.value.value }.toString()
        }

    }

    companion object {
        fun identity(size: Int): HashRowColIntMatrix {
            val I = HashRowColIntMatrix(size,size)
            for(i in 0 until size) {
                I[i,i] = 1
            }
            return I
        }
    }


}