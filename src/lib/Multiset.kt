package lib

import java.io.Serializable
import kotlin.math.max
import kotlin.math.min

interface Multiset<T>: Set<T> {
    val counts: Map<T,Int>

    val supportSet: Set<T>
        get() = counts.keys

    fun count(m : T) : Int {
        return counts.getOrDefault(m,0)
    }

    fun union(other : Multiset<T>) : MutableMultiset<T> {
        val result = other.toMutableMultiset()
        this.counts.forEach { result.merge(it.key, it.value, ::max) }
        return result
    }

    fun union(other : Iterable<T>) : MutableMultiset<T> {
        val result = other.toMutableMultiset()
        this.counts.forEach { result.merge(it.key, it.value, ::max) }
        return result
    }

    fun intersect(other : Iterable<T>) : MutableMultiset<T> {
        val result = HashMultiset<T>()
        other.forEach {
            if(count(it) > result.count(it)) result.add(it)
        }
        return result
    }

    fun intersect(other : Multiset<T>) : MutableMultiset<T> {
        val result = HashMultiset<T>()
        other.counts.forEach {
            val intersect = min(this.count(it.key), it.value)
            if(intersect > 0) result.add(it.key, intersect)
        }
        return result
    }


    operator fun minus(other: Multiset<T>): MutableMultiset<T> {
        val result = HashMultiset<T>()
        counts.forEach { entry ->
            val newCount = entry.value - other.count(entry.key)
            if(newCount > 0) result.add(entry.key, newCount)
        }
        return result
    }

    operator fun minus(other: Iterable<T>) = subtract(other)

    fun subtract(other: Iterable<T>): MutableMultiset<T> {
        val result = this.toMutableMultiset()
        other.forEach { result.remove(it) }
        return result
    }


    operator fun plus(other: Iterable<T>): MutableMultiset<T> {
        val result = this.toMutableMultiset()
        other.forEach { result.add(it) }
        return result
    }

    override fun contains(element: T): Boolean {
        return counts.containsKey(element)
    }


    fun isSubsetOf(other: Multiset<T>): Boolean {
        this.counts.entries.forEach { if (it.value > other.count(it.key)) return false }
        return true
    }


    override fun containsAll(elements: Collection<T>): Boolean {
        return counts.keys.containsAll(elements)
    }

    override fun isEmpty(): Boolean {
        return counts.isEmpty()
    }


}

fun <T>emptyMultiset(): Multiset<T> = emptySet<T>().asMultiSet()

fun <T>Set<T>.asMultiSet(): Multiset<T> {
    val set = this
    return object: Serializable, Multiset<T> {
        override val counts: Map<T,Int>
            get() = set.associateWith { 1 }
        override val size: Int
            get() = set.size

        override fun count(m: T)= if(set.contains(m)) 1 else 0
        override fun contains(element: T) = set.contains(element)
        override fun containsAll(elements: Collection<T>) = set.containsAll(elements)
        override fun isEmpty() = set.isEmpty()
        override fun iterator() = set.iterator()
    }
}


fun <T>Collection<T>.toMultiset(): Multiset<T> {
    return this.toMutableMultiset()
}

fun <T>Iterable<T>.toMultiset(): Multiset<T> {
    return this.toMutableMultiset()
}

fun <T>Sequence<T>.toMultiset(): Multiset<T> {
    return this.toMutableMultiset()
}

fun <T>Array<out T>.toMultiset(): Multiset<T> {
    return this.toMutableMultiset()
}
