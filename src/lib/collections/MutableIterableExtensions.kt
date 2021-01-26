package lib.collections


inline fun<T> MutableIterable<T>.mutableWithIndex(): MutableIterable<IndexedValue<T>> = MutableIterable {
    IndexingMutableIterator(this)
}

inline fun<T> Iterable<MutableIterable<T>>.mutableFlatten(): MutableIterable<T> = MutableIterable {
    FlattenedMutableIterator(this)
}

fun<T,R> Iterable<T>.mutableFlatMap(transform: (T)->MutableIterable<R>): MutableIterable<R> = MutableIterable {
    MutableFlatMapIterator(this, transform)
}

fun<T,R> MutableIterable<T>.mutableMap(transform: (T)->R): MutableIterable<R> = MutableIterable {
    MutableMapIterator(this, transform)
}

fun<T,R> MutableIterable<T>.mutableMapByRemovable(transform: (RemovableEntry<T>)->R): MutableIterable<R> = MutableIterable {
    MutableMapByRemovableIterator(this, transform)
}


inline fun<T> MutableIterable(crossinline iteratorFactory: ()->MutableIterator<T>): MutableIterable<T> {
    return object: MutableIterable<T> {
        override fun iterator() = iteratorFactory()
    }
}

class FlattenedMutableIterator<T>: MutableIterator<T> {
    val outerIterator: Iterator<MutableIterable<T>>
    var innerIterator: MutableIterator<T>
    var removeIterator: MutableIterator<T>? = null

    constructor(iterableOfIterables: Iterable<MutableIterable<T>>) {
        outerIterator = iterableOfIterables.iterator()
        if(outerIterator.hasNext()) {
            innerIterator = outerIterator.next().iterator()
        } else {
            innerIterator = mutableListOf<T>().iterator()
        }
    }

    override fun hasNext() = moveToNextElement()

    override fun next(): T {
        moveToNextElement()
        removeIterator = innerIterator
        return innerIterator.next()
    }

    override fun remove() { removeIterator?.remove()?:throw(NoSuchElementException()) }

    private fun moveToNextElement(): Boolean {
        if(innerIterator.hasNext()) return true
        while(outerIterator.hasNext()) {
            innerIterator = outerIterator.next().iterator()
            if(innerIterator.hasNext()) return true
        }
        return false
    }
}


class MutableFlatMapIterator<T,R>: MutableIterator<R> {
    val outerIterator: Iterator<T>
    val elementToIterable: (T)->MutableIterable<R>
    var innerIterator: MutableIterator<R>
    var removeIterator: MutableIterator<R>? = null

    constructor(iterableOfIterables: Iterable<T>, elementToIterable: (T)->MutableIterable<R>) {
        this.elementToIterable = elementToIterable
        outerIterator = iterableOfIterables.iterator()
        if(outerIterator.hasNext()) {
            innerIterator = elementToIterable(outerIterator.next()).iterator()
        } else {
            innerIterator = mutableListOf<R>().iterator()
        }
    }

    override fun hasNext() = moveToNextElement()

    override fun next(): R {
        moveToNextElement()
        removeIterator = innerIterator
        return innerIterator.next()
    }

    override fun remove() { removeIterator?.remove()?:throw(NoSuchElementException()) }

    private fun moveToNextElement(): Boolean {
        if(innerIterator.hasNext()) return true
        while(outerIterator.hasNext()) {
            innerIterator = elementToIterable(outerIterator.next()).iterator()
            if(innerIterator.hasNext()) return true
        }
        return false
    }
}

class MutableMapIterator<T,R>(val teeIterator: MutableIterator<T>, val transform: (T)->R) : MutableIterator<R> {
    constructor(teeIterable: MutableIterable<T>, transform: (T)->R): this(teeIterable.iterator(), transform)
    override fun hasNext() = teeIterator.hasNext()
    override fun next()= transform(teeIterator.next())
    override fun remove() { teeIterator.remove() }
}

class MutableMapByRemovableIterator<T,R>(val teeIterator: MutableIterator<T>, val transform: (RemovableEntry<T>)->R) : MutableIterator<R> {
    constructor(teeIterable: MutableIterable<T>, transform: (RemovableEntry<T>)->R): this(teeIterable.iterator(), transform)
    override fun hasNext() = teeIterator.hasNext()
    override fun next()= transform(RemovableEntry(teeIterator.next(),teeIterator::remove))
    override fun remove() { teeIterator.remove() }
}


class IndexingMutableIterator<out T>(private val iterator: MutableIterator<T>) : MutableIterator<IndexedValue<T>> {
    constructor(iterable: MutableIterable<T>): this(iterable.iterator())
    private var index = 0
    override fun hasNext(): Boolean = iterator.hasNext()
    override fun next(): IndexedValue<T> = IndexedValue(index++, iterator.next())
    override fun remove() {
        iterator.remove()
    }
}

class RemovableEntry<T>(val value: T, val remove: ()->Unit)
