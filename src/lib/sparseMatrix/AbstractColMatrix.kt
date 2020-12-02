package lib.sparseMatrix

import lib.vector.MutableSparseVector
import lib.vector.SparseVector

abstract class AbstractColMatrix<T: Any>(override val nRows: Int, override val columns: MutableList<MutableSparseVector<T>>): SparseColMatrix<T> {

    fun resize(nRows: Int, nCols: Int) {
        if(nCols > columns.size) {
            repeat(nCols - columns.size) { columns.add(columns.first().new()) }
        } else if(nCols < columns.size) {
            repeat(columns.size - nCols) {
                columns.removeAt(columns.lastIndex)
            }
        }
        if(nRows < this.nRows) {
            for(j in columns.indices) {
                val it = columns[j].nonZeroEntries.iterator()
                while(it.hasNext()) {
                    val entry = it.next()
                    if(entry.key >= nRows) it.remove()
                }
            }
        }
    }

}