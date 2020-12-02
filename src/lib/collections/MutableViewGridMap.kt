package lib.collections

import java.lang.RuntimeException
import java.lang.StringBuilder
import java.util.function.BiFunction

class MutableViewGridMap<T>(
    val _columnData: ArrayList<MutableMap<Int,T>>,
    val _rowData: ArrayList<MutableSet<Int>>
    ) {

    inline val nCols: Int
        get() = _columnData.size

    inline val nRows: Int
        get() = _rowData.size

    val columns: List<MutableMap<Int,T>> = ColListView()

    val rows: List<MutableMap<Int,T>> = RowListView()


    constructor(): this(ArrayList<MutableMap<Int,T>>(), ArrayList<MutableSet<Int>>())

    constructor(nRows: Int, nCols: Int): this(ArrayList<MutableMap<Int,T>>(nCols), ArrayList<MutableSet<Int>>(nRows)) {
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


    inline fun compute(row: Int, col: Int, crossinline remappingFunction: (Int, T?)->T?): T? {
        return _columnData[col].compute(row) { row, oldValue ->
            val newValue = remappingFunction(row, oldValue)
            if(oldValue == null) {
                if(newValue != null) _rowData[row].add(col)
            } else if (newValue == null) _rowData[row].remove(col)
            newValue
        }
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

//    inner class RowListView: AbstractList<Map<Int,T>>() {
//        override val size: Int
//            get() = _rowData.size
//
//        override fun get(index: Int): Map<Int, T> {
//            return RowMapView(index)
//        }
//    }

//    inner class EntrySet(val row: Int): AbstractSet<Map.Entry<Int, T>>() {
//        override val size: Int
//            get() = _rowData[row].size
//
//        override fun iterator() = IteratorWrapper(row, _rowData[row].iterator())
//    }
//
//
//    inner class IteratorWrapper(val row: Int, val colIterator: Iterator<Int>) : Iterator<Map.Entry<Int, T>> {
//        override fun hasNext() = colIterator.hasNext()
//
//        override fun next(): Map.Entry<Int, T> {
//            val colIndex = colIterator.next()
//            return java.util.AbstractMap.SimpleEntry<Int,T>(colIndex, _columnData[colIndex][row])
//        }
//    }

    inner class RowListView: AbstractList<MutableMap<Int,T>>() {
        override val size: Int
            get() = _rowData.size
        override fun get(index: Int) = RowMapView(index)
    }


    inner class ColListView: AbstractList<MutableMap<Int,T>>() {
        override val size: Int
            get() = _columnData.size
        override fun get(index: Int) = ColMapView(index)
    }


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
            return this@MutableViewGridMap.compute(key, col, remappingFunction)
        }
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
            return this@MutableViewGridMap.compute(row, key, remappingFunction)
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
            _columnData[lastSeenCol?:throw(RuntimeException("Can't remove unseen item"))].remove(row)
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