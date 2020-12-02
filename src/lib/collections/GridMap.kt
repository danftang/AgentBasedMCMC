package lib.collections

import java.lang.StringBuilder

class GridMap<T>(val _columnData: ArrayList<MutableMap<Int,T>>, val _rowData: ArrayList<MutableSet<Int>>) {

    inline val nCols: Int
        get() = _columnData.size

    inline val nRows: Int
        get() = _rowData.size

    inline val columns: List<Map<Int,T>>
        get() = _columnData

    val rows: List<Map<Int,T>> = RowListView()

    constructor(): this(ArrayList<MutableMap<Int,T>>(), ArrayList<MutableSet<Int>>())

    constructor(nRows: Int, nCols: Int): this(ArrayList<MutableMap<Int,T>>(nCols), ArrayList<MutableSet<Int>>(nRows)) {
        repeat(nCols) { addColumn() }
        repeat(nRows) { addRow() }
    }

    operator fun get(row: Int, col: Int): T? {
        return _columnData[col][row]
    }

    operator fun set(row: Int, col: Int, value: T) {
        _columnData[col].compute(row) { row, oldValue ->
            if(oldValue == null) _rowData[row].add(col)
            value
        }
    }

    fun remove(row: Int, col: Int) {
        _columnData[col].compute(row) { row, oldValue ->
            if(oldValue != null) _rowData[row].remove(col)
            null
        }
    }

//    fun computeNonNullMapping(row: Int, col: Int, remappingFunction: (Int, T?)->T) {
//        _columnData[col].compute(row, remappingFunction)
//    }

    fun resize(nRows: Int, nCols: Int) {
        if(nRows > this.nRows) {
            repeat(nRows - this.nRows) { addRow() }
        } else if(nRows < this.nRows) {
            repeat(this.nRows - nRows) {
                val rowIndex = _rowData.lastIndex
                _rowData[rowIndex].forEach { colIndex ->
                    _columnData[colIndex].remove(rowIndex)
                }
                _rowData.removeAt(rowIndex)
            }
        }
        if(nCols > this.nCols) {
            repeat(nCols - this.nCols) { addColumn() }
        } else if(nCols < this.nCols) {
            repeat(this.nCols - nCols) {
                val colIndex = _columnData.lastIndex
                _columnData[colIndex].keys.forEach { rowIndex ->
                    _rowData[rowIndex].remove(colIndex)
                }
                _columnData.removeAt(colIndex)
            }
        }
    }

    fun addRow() {
        _rowData.add(HashSet(4))
    }

    fun addColumn() {
        _columnData.add(HashMap(4))
    }

    fun columnReassign(col: Int, remappingFunction: (T) -> T) {
        val colMap = _columnData[col]
        for(entry in colMap) {
            entry.setValue(remappingFunction(entry.value))
        }
    }


    inline fun rowReassign(row: Int, crossinline remappingFunction: (T) -> T) {
        val rowSet = _rowData[row]
        for(col in rowSet) {
            _columnData[col].compute(row) { _, oldValue ->
                remappingFunction(oldValue!!)
            }
        }
    }


    inline fun compute(row: Int, col: Int, crossinline remappingFunction: (Int, T?)->T?) {
        _columnData[col].compute(row) { row, oldValue ->
            val newValue = remappingFunction(row, oldValue)
            if(oldValue == null) {
                if(newValue != null) _rowData[row].add(col)
            } else {
                if (newValue == null) _rowData[row].remove(col)
            }
            newValue
        }
    }

    fun merge(row: Int, col: Int, value: T, remappingFunction: (T,T)->T?) {
        _columnData[col].compute(row) { row, oldValue ->
            if(oldValue == null) {
                _rowData[row].add(col)
                value
            } else {
                val mappedValue = remappingFunction(oldValue, value)
                if(mappedValue == null) _rowData[row].remove(col)
                mappedValue
            }
        }
    }

    override fun toString(): String {
        val out = StringBuilder()
        for(row in 0 until nRows) {
            for(col in 0 until nCols) {
                val v = get(row,col)
                val vs = v.toString().padEnd(6,' ').take(6)
                if(v == null) out.append("......\t") else out.append("$vs\t")
            }
            out.append('\n')
        }
        return out.toString()
    }

    inner class RowListView: AbstractList<Map<Int,T>>() {
        override val size: Int
            get() = _rowData.size

        override fun get(index: Int): Map<Int, T> {
            return RowMapView(index)
        }
    }


    inner class RowMapView(val row: Int): AbstractMap<Int,T>() {

        override val entries: Set<Map.Entry<Int, T>> = EntrySet(row)

        override fun get(col: Int): T? {
            return _columnData[col][row]
        }

        override fun containsKey(col: Int): Boolean {
            return _rowData[row].contains(col)
        }
    }


    inner class EntrySet(val row: Int): AbstractSet<Map.Entry<Int, T>>() {
        override val size: Int
            get() = _rowData[row].size

        override fun iterator() = IteratorWrapper(row, _rowData[row].iterator())
    }


    inner class IteratorWrapper(val row: Int, val colIterator: Iterator<Int>) : Iterator<Map.Entry<Int, T>> {
        override fun hasNext() = colIterator.hasNext()

        override fun next(): Map.Entry<Int, T> {
            val colIndex = colIterator.next()
            return java.util.AbstractMap.SimpleEntry<Int,T>(colIndex, _columnData[colIndex][row])
        }
    }

}