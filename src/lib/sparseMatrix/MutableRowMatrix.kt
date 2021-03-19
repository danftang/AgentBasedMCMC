package lib.sparseMatrix

import lib.collections.RemovableEntry
import lib.collections.mutableMapByRemovable
import lib.sparseVector.MutableSparseVector
import lib.sparseVector.SparseVector

interface MutableRowMatrix<T>: RowMatrix<T>, MutableMatrix<T> {


//    override val rows: List<MutableSparseVector<T>>
//    override val nonZeroEntries: Iterable<EntryMatrix.Entry<T>>
//        get() = rows.withIndex().flatMap { (row, vector) ->
//            vector.nonZeroEntries.entries.mutableMapByRemovable { entry ->
//                REntry(row, entry) { this.isZero() }
//            }
//        }
//

//    override fun mapAssign(row: Int, col: Int, remappingFunction: (T) -> T) {
//        rows[row].mapAssign(col, remappingFunction)
//    }
//
//
//
//    class REntry<T>(override val row: Int, val entry: RemovableEntry<MutableMap.MutableEntry<Int, T>>, val isZero: T.()->Boolean): EntryMatrix.Entry<T> {
//        override val col: Int
//            get() = entry.value.key
//        override val value: T
//            get() = entry.value.value
//        override fun setValue(newValue: T): T {
//            return if(newValue.isZero()) {
//                val oldValue = value
//                entry.remove()
//                oldValue
//            } else {
//                entry.value.setValue(newValue)
//            }
//        }
//    }
}