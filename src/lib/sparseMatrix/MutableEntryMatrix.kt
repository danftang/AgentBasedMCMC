package lib.sparseMatrix

interface MutableEntryMatrix<T>: EntryMatrix<T>, MutableMatrix<T> {
//    interface MutableEntry<T>: EntryMatrix.Entry<T> {
//        fun setValue(newValue: T): T
//    }

//    override val nonZeroEntries: MutableIterable<MutableEntry<T>>

}