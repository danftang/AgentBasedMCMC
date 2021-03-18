package lib.sparseMatrix

interface MutableEntryMatrix<T>: EntryMatrix<T> {
    interface MutableEntry<T>: EntryMatrix.Entry<T> {
        fun setValue(newValue: T): T
    }

    override val nonZeroEntries: MutableIterable<MutableEntry<T>>
    operator fun set(row: Int, col: Int, value: T)
}