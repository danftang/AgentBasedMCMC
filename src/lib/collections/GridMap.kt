package lib.collections

import java.lang.RuntimeException
import java.lang.StringBuilder
import java.util.AbstractMap

class GridMap<T>(
    private val _columnData: ArrayList<MutableMap<Int,T>>,
    private val _rowData: ArrayList<MutableSet<Int>>
    ) {

    private val columnAccessors = ArrayList<ColMapView>()
    private val rowAccessors    = ArrayList<RowMapView>()

    val nCols: Int
        get() = _columnData.size

    val nRows: Int
        get() = _rowData.size

    val size: Int
        get() = _rowData.sumBy { it.size }

    val columns: List<MutableMap<Int,T>>
        get() = columnAccessors

    val rows: List<MutableMap<Int,T>>
        get() = rowAccessors

    val entries: Sequence<GridMapEntry<T>>
        get() = _columnData.asSequence().withIndex().flatMap { column ->
            column.value.asSequence().map { colEntry ->
                GridMapEntry(colEntry.key, column.index, colEntry.value)
            }
        }

    data class GridMapEntry<T>(val row: Int, val col: Int, val value: T)

    constructor(nRows: Int, nCols: Int): this(ArrayList(nCols), ArrayList(nRows)) {
        repeat(nCols) { addColumn() }
        repeat(nRows) { addRow() }
    }

    operator fun get(row: Int, col: Int): T? {
        return _columnData[col][row]
    }

    operator fun set(row: Int, col: Int, value: T) {
        if(_columnData[col].put(row, value) == null) _rowData[row].add(col)
    }

    fun remove(row: Int, col: Int) {
        if(_columnData[col].remove(row) != null) _rowData[row].remove(col)
    }

    fun clear() {
        for(colMap in _columnData) colMap.clear()
        for(rowMap in _rowData) rowMap.clear()
    }


//    fun computeNonNullMapping(row: Int, col: Int, remappingFunction: (Int, T?)->T) {
//        _columnData[col].compute(row, remappingFunction)
//    }

    fun resize(nRows: Int, nCols: Int) {
        if(nRows > this.nRows) {
            repeat(nRows - this.nRows) { addRow() }
        } else if(nRows < this.nRows) {
            repeat(this.nRows - nRows) { removeRow() }
        }
        if(nCols > this.nCols) {
            repeat(nCols - this.nCols) { addColumn() }
        } else if(nCols < this.nCols) {
            repeat(this.nCols - nCols) { removeColumn() }
        }
    }

    fun addRow() {
        _rowData.add(HashSet())
        rowAccessors.add(RowMapView(_rowData.lastIndex))
    }

    fun addColumn() {
        _columnData.add(HashMap())
        columnAccessors.add(ColMapView(_columnData.lastIndex))
    }

    fun removeRow() {
        val rowIndex = _rowData.lastIndex
        _rowData[rowIndex].forEach { colIndex ->
            _columnData[colIndex].remove(rowIndex)
        }
        _rowData.removeAt(rowIndex)
        rowAccessors.removeAt(rowIndex)
    }

    fun removeColumn() {
        val colIndex = _columnData.lastIndex
        _columnData[colIndex].keys.forEach { rowIndex ->
            _rowData[rowIndex].remove(colIndex)
        }
        _columnData.removeAt(colIndex)
        columnAccessors.removeAt(colIndex)
    }

    fun columnReassign(col: Int, remappingFunction: (T) -> T) {
        val colMap = _columnData[col]
        for(entry in colMap) {
            entry.setValue(remappingFunction(entry.value))
        }
    }


    fun rowReassign(row: Int, remappingFunction: (T) -> T) {
        val rowSet = _rowData[row]
        for(col in rowSet) {
            _columnData[col].compute(row) { _, oldValue ->
                remappingFunction(oldValue!!)
            }
        }
    }


    fun compute(row: Int, col: Int, remappingFunction: (Int, T?)->T?): T? {
        return _columnData[col].compute(row) { rowIndex, oldValue ->
            val newValue = remappingFunction(rowIndex, oldValue)
            if(oldValue == null) {
                if(newValue != null) _rowData[rowIndex].add(col)
            } else if (newValue == null) _rowData[rowIndex].remove(col)
            newValue
        }
    }

    fun swapRows(row1: Int, row2: Int) {
        for(entry in rows[row1]) {
            columns[entry.key][row2] = entry.value
            columns[entry.key].remove(row1)
        }
        for(entry in rows[row2]) {
            columns[entry.key][row1] = entry.value
            if(!rows[row1].containsKey(entry.key)) columns[entry.key].remove(row2)
        }
        val tmpRow = _rowData[row1]
        _rowData[row1] = _rowData[row2]
        _rowData[row2] = tmpRow
    }


//    fun merge(row: Int, col: Int, value: T, remappingFunction: (T,T)->T?) {
//        _columnData[col].compute(row) { row, oldValue ->
//            if(oldValue == null) {
//                _rowData[row].add(col)
//                value
//            } else {
//                val mappedValue = remappingFunction(oldValue, value)
//                if(mappedValue == null) _rowData[row].remove(col)
//                mappedValue
//            }
//        }
//    }

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


//    inner class RowListView: AbstractList<MutableMap<Int,T>>() {
//        override val size: Int
//            get() = _rowData.size
//        override fun get(index: Int) = RowMapView(index)
//    }
//
//
//    inner class ColListView: AbstractList<MutableMap<Int,T>>() {
//        override val size: Int
//            get() = _columnData.size
//        override fun get(index: Int) = ColMapView(index)
//    }


    inner class ColMapView(val col: Int): MutableMap<Int, T> by _columnData[col] {
        override val entries: MutableSet<MutableMap.MutableEntry<Int, T>> = ColEntrySet(col)

        override fun remove(key: Int): T? {
            val oldVal = _columnData[col].remove(key)
            if(oldVal != null) _rowData[key].remove(col)
            return oldVal
        }

        override fun put(key: Int, value: T): T? {
            val oldValue = _columnData[col].put(key,value)
            if(oldValue == null) _rowData[key].add(col)
            return oldValue
        }

        fun compute(key: Int, remappingFunction: (Int, T?)->T): T? {
            return this@GridMap.compute(key, col, remappingFunction)
        }

        override fun toString() = entries.toString()
    }

    inner class RowMapView(val row: Int): AbstractMutableMap<Int,T>() {
        override val entries: MutableSet<MutableMap.MutableEntry<Int, T>> = RowEntrySet(row)

        override fun get(key: Int): T? {
            return _columnData[key][row]
        }

        override fun remove(key: Int):T? {
            val oldVal = _columnData[key].remove(row)
            if(oldVal != null) _rowData[row].remove(key)
            return oldVal
        }

        override fun containsKey(key: Int): Boolean {
            return _rowData[row].contains(key)
        }

        override fun put(key: Int, value: T): T? {
            val oldValue = _columnData[key].put(row, value)
            if(oldValue == null) _rowData[row].add(key)
            return oldValue
        }

        fun compute(key: Int, remappingFunction: (Int, T?)->T): T? {
            return this@GridMap.compute(row, key, remappingFunction)
        }


    }


    inner class ColEntrySet(val col: Int): AbstractMutableSet<MutableMap.MutableEntry<Int, T>>() {
        override val size: Int
            get() = _columnData[col].size

        override fun add(element: MutableMap.MutableEntry<Int, T>): Boolean {
            _columnData[col][element.key] = element.value
            _rowData[element.key].add(col)
            return true
        }

        override fun iterator() = ColIteratorWrapper(col, _columnData[col].iterator())

        override fun toString(): String {
            val out = StringBuilder()
            val it = iterator()
            out.append('[')
            while(it.hasNext()) {
                out.append(it.next().toString())
                out.append(',')
            }
            out.deleteCharAt(out.lastIndex)
            out.append(']')
            return out.toString()
        }
    }

    inner class RowEntrySet(val row: Int): AbstractMutableSet<MutableMap.MutableEntry<Int, T>>() {
        override val size: Int
            get() = _rowData[row].size

        override fun add(element: MutableMap.MutableEntry<Int, T>): Boolean {
            _columnData[element.key][row] = element.value
            _rowData[row].add(element.key)
            return true
        }

        override fun iterator() = RowIteratorWrapper(row, _rowData[row].iterator())
    }

    inner class ColIteratorWrapper(val col: Int, val underlyingIterator: MutableIterator<MutableMap.MutableEntry<Int, T>>)
        : MutableIterator<MutableMap.MutableEntry<Int, T>> {
        var lastSeenRow: Int? = null

        override fun hasNext() = underlyingIterator.hasNext()

        override fun next(): MutableMap.MutableEntry<Int, T> {
            val underlyingEntry = underlyingIterator.next()
            lastSeenRow = underlyingEntry.key
//            assert(lastSeenRow != null)
            return underlyingEntry
        }

        override fun remove() {
            _rowData[lastSeenRow?:throw(RuntimeException("Can't remove unseen item"))].remove(col)
            underlyingIterator.remove()
        }
    }

    inner class RowIteratorWrapper(val row: Int, val underlyingIterator: MutableIterator<Int>)
        : MutableIterator<MutableMap.MutableEntry<Int, T>> {
        var lastSeenCol: Int? = null
        override fun hasNext() = underlyingIterator.hasNext()

        override fun next(): MutableMap.MutableEntry<Int, T> {
            val col = underlyingIterator.next()
            lastSeenCol = col
            return RowEntry(row, col)
        }

        override fun remove() {
            _columnData[lastSeenCol?:throw(RuntimeException("Trying to remove from an iterator before calling next()"))].remove(row)
            underlyingIterator.remove()
        }
    }


    inner class RowEntry(val row: Int, val col: Int): MutableMap.MutableEntry<Int, T> {
        override val key: Int
            get() = col
        override val value: T
            get() = _columnData[col][row]!!

        override fun setValue(newValue: T): T {
            return _columnData[col].put(row, newValue)!!
        }
    }
}