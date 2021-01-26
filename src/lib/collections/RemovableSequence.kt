package lib.collections

import java.lang.RuntimeException

fun<T> MutableIterable<T>.asSequenceOfRemovables(): Sequence<Removable<T>> = MutableSequence(this)
fun<T> Sequence<Removable<T>>.asMutableIterable(): MutableIterable<T> = RemovableToMutableIterable(this)

class Removable<T>(val value: T, val remove: ()->Unit)

inline class MutableToRemovableIterator<T>(val it: MutableIterator<T>): Iterator<Removable<T>> {
    override fun hasNext() = it.hasNext()
    override fun next() = Removable(it.next(), it::remove)
}

inline class MutableSequence<T>(val data: MutableIterable<T>): Sequence<Removable<T>> {
    override fun iterator(): Iterator<Removable<T>> = MutableToRemovableIterator(data.iterator())
}

class RemovableToMutableIterator<T>(val it: Iterator<Removable<T>>): MutableIterator<T> {
    var last: Removable<T>? = null

    override fun hasNext() = it.hasNext()
    override fun next(): T {
        last = it.next()
        return last!!.value
    }

    override fun remove() {
        last?.remove?.invoke()?:throw(RuntimeException("Can't call remove() before calling next()"))
    }

}

inline class RemovableToMutableIterable<T>(val mSequenece: Sequence<Removable<T>>): MutableIterable<T> {
    override fun iterator(): MutableIterator<T> = RemovableToMutableIterator(mSequenece.iterator())
}