package lib.collections

import java.lang.RuntimeException

// Creates a mutable view of a mutable map, where the values are transformed.
// override the removeHook function if you want to do any processing when
// an element is removed from the map.
abstract class MutableMapView<K,UNDERLYINGV, V>(val underlyingMap: MutableMap<K,UNDERLYINGV>): AbstractMutableMap<K,V>() {

    override val entries: MutableSet<MutableMap.MutableEntry<K, V>> = EntrySet()

    inner class EntrySet: AbstractMutableSet<MutableMap.MutableEntry<K, V>>() {
        override val size: Int
            get() = underlyingMap.size

        override fun add(element: MutableMap.MutableEntry<K, V>): Boolean {
            put(element.key, element.value)
            return true
        }

        override fun iterator() = IteratorWrapper(underlyingMap.iterator())
    }


    inner class IteratorWrapper(
        val underlyingIterator: MutableIterator<MutableMap.MutableEntry<K, UNDERLYINGV>>,
        var lastSeenKey: K? = null
    ) : MutableIterator<MutableMap.MutableEntry<K, V>> {
        override fun hasNext() = underlyingIterator.hasNext()

        override fun next(): MutableMap.MutableEntry<K, V> {
            val underlyingEntry = underlyingIterator.next()
            lastSeenKey = underlyingEntry.key
            return Entry(underlyingEntry)
        }

        override fun remove() {
            this@MutableMapView.removeHook(lastSeenKey?:throw(RuntimeException("No element to remove")))
            underlyingIterator.remove()
        }

    }


    inner class Entry(val underlyingEntry: MutableMap.MutableEntry<K,UNDERLYINGV>): MutableMap.MutableEntry<K, V> {
        override val key: K
        get() = underlyingEntry.key
        override val value: V
        get() = valueTransform(underlyingEntry.value)

        override fun setValue(newValue: V): V {
            set(underlyingEntry, newValue)
            return newValue
        }

    }

//    override fun put(key: K, value: V): V?

    // Sets a mutable entry to the given newValue
    abstract fun set(entry: MutableMap.MutableEntry<K, UNDERLYINGV>, newValue: V)

    // Transforms a value in the underlying map to a value in the view
    abstract fun valueTransform(underlyingValue: UNDERLYINGV): V

    // override this function to perform any processing on removal of an element
    open fun removeHook(key: K) {
    }


    override fun get(key: K): V? {
        return underlyingMap[key]?.let { valueTransform(it) }
    }

    override fun remove(key: K):V? {
        removeHook(key)
        return underlyingMap.remove(key)?.let { valueTransform(it) }
    }

    override fun containsKey(key: K): Boolean {
        return underlyingMap.containsKey(key)
    }
}