package lib

import java.io.Serializable

open class HashMultiset<T>(private val mutableCounts: HashMap<T,Int> = HashMap()): MutableMultiset<T>, Serializable {
    override val counts: Map<T,Int>
        get() = mutableCounts
    override var size : Int = 0


    constructor(initialCapacity : Int) : this(HashMap(initialCapacity))

    constructor(container : Sequence<T>) : this() {
        container.forEach { add(it) }
    }

    constructor(container : Iterable<T>) : this() {
        container.forEach { add(it) }
    }

    constructor(container : Collection<T>) : this(HashMap(container.size)) {
        container.forEach { add(it) }
    }

    constructor(array: Array<out T>): this(array.asList())

    override fun compute(member: T, transform: (T, Int) -> Int): Int {
        val newVal = mutableCounts.compute(member) { member, oldValOrNull ->
            val oldVal = oldValOrNull?:0
            val newVal = transform(member, oldVal)
            if(newVal <= 0) {
                size -= oldVal
                null
            } else {
                size += newVal - oldVal
                newVal
            }
        }
        return newVal?:0
    }


    override fun merge(member: T, value: Int, transform: (Int, Int) -> Int): Int {
        size += value
        val newVal = mutableCounts.merge(member, value) { oldVal, otherVal ->
            val newVal = transform(oldVal,otherVal)
            if(newVal <= 0) {
                size -= oldVal + otherVal
                null
            } else {
                size += newVal - oldVal - otherVal
                newVal
            }
        }
        return newVal?:0
    }




    override fun retainAll(elements: Multiset<T>): Boolean {
        return mutableCounts.entries.removeIf { entry ->
            val toRetain = elements.count(entry.key)
            if (toRetain < entry.value) {
                size -= entry.value
                if (toRetain > 0) {
                    entry.setValue(toRetain)
                    size += toRetain
                    false
                } else {
                    true
                }
            }
            false
        }
    }



    override operator fun iterator() : MutableIterator<T> {
        return MultiMutableIterator(mutableCounts.iterator(), this::decrementSize)
    }




//    override fun equals(other: Any?): Boolean {
//        if(other is lib.HashMultiset<*>) {
//            return this.map == other.map
//        }
//        return false
//    }


    private fun decrementSize() {--size}

    ////////////////////////////////////////////////////////
    // overrides of standard algorithms that are quicker
    // than the default implementations
    ////////////////////////////////////////////////////////
    fun elementAt(index : Int) : T {
        var count = index
        val it = mutableCounts.iterator()
        var entry : MutableMap.MutableEntry<T,Int>
        do {
            entry = it.next()
            count -= entry.value
        } while(count>=0)
        return entry.key
    }


    override fun equals(other: Any?): Boolean {
        if(this === other) return true
        if(other is Multiset<*>) {
            return this.counts == other.counts
        }
        if(other is Set<*>) {
            if(!this.counts.filter { it.value != 1 }.isEmpty()) return false
            return this.supportSet == other
        }
        return false
    }




    class MultiMutableIterator<A>(val mapIterator : MutableIterator<MutableMap.MutableEntry<A,Int>>, val decrementSize : ()->Unit) : MutableIterator<A> {
        private var currentEntry : MutableMap.MutableEntry<A,Int>? = null
        private var nItems : Int = 0

        override fun remove() {
            val entry = currentEntry
            when {
                entry == null       -> throw(NoSuchElementException())
                entry.value == 1    -> mapIterator.remove()
                else                -> entry.setValue(entry.value - 1)
            }
            decrementSize()
        }

        override fun hasNext(): Boolean = mapIterator.hasNext() || nItems > 0

        override fun next(): A {
            if(nItems > 0) {
                --nItems
            } else {
                val nextEntry = mapIterator.next()
                currentEntry = nextEntry
                nItems = nextEntry.value-1
            }
            return currentEntry!!.key
        }

    }

    override fun hashCode(): Int {
        return mutableCounts.hashCode()
    }

    override fun toString(): String {
        //  String conversion
        val it: Iterator<T> = this.iterator()
        if (this.isEmpty()) return "[]"
        val sb = StringBuilder()
        sb.append('[')
        while (true) {
            val e = it.next()
            sb.append(if (e === this) "(this Collection)" else e)
            if (!it.hasNext()) return sb.append(']').toString()
            sb.append(',').append(' ')
        }
    }

    override fun clear() {
        mutableCounts.clear()
    }


}


fun <T>hashMultisetOf(vararg elements : T) : HashMultiset<T> {
    return HashMultiset(elements.asList())
}

