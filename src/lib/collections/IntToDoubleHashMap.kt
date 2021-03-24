package lib.collections

import org.apache.commons.math3.util.OpenIntToDoubleHashMap

class IntToDoubleHashMap(val apacheMap: OpenIntToDoubleHashMap): AbstractMutableMap<Int,Double>() {
    override val size: Int
        get() = apacheMap.size()

    override fun containsKey(key: Int) = apacheMap.containsKey(key)

    override fun get(key: Int): Double {
        return apacheMap.get(key)
    }

    operator fun set(key: Int, value: Double) {
        apacheMap.put(key, value)
    }

    override fun put(key: Int, value: Double): Double {
        return apacheMap.put(key, value)
    }
//
//    override fun isEmpty(): Boolean {
//        return apacheMap.size() == 0
//    }

    override val entries: MutableSet<MutableMap.MutableEntry<Int, Double>> = EntrySet(apacheMap)
//    override val keys: MutableSet<Int>
//        get() = TODO("Not yet implemented")
//    override val values: MutableCollection<Double>
//        get() = TODO("Not yet implemented")

//    override fun clear() {
//        TODO("Not yet implemented")
//    }
//
//    override fun containsValue(value: Double): Boolean {
//        TODO("Not yet implemented")
//    }
//
//
//    override fun putAll(from: Map<out Int, Double>) {
//        TODO("Not yet implemented")
//    }
//
//    override fun remove(key: Int): Double {
//        return apacheMap.remove(key)
//    }

    inline class EntrySet(val apacheMap: OpenIntToDoubleHashMap): MutableSet<MutableMap.MutableEntry<Int,Double>> {
        override val size: Int
            get() = apacheMap.size()

        override fun add(element: MutableMap.MutableEntry<Int, Double>): Boolean {
            TODO("Not yet implemented")
        }

        override fun iterator(): MutableIterator<MutableMap.MutableEntry<Int, Double>> {
            return IteratorWrapper(apacheMap.iterator())
        }

        override fun addAll(elements: Collection<MutableMap.MutableEntry<Int, Double>>): Boolean {
            TODO("Not yet implemented")
        }

        override fun clear() {
            TODO("Not yet implemented")
        }

        override fun remove(element: MutableMap.MutableEntry<Int, Double>): Boolean {
            TODO("Not yet implemented")
        }

        override fun removeAll(elements: Collection<MutableMap.MutableEntry<Int, Double>>): Boolean {
            TODO("Not yet implemented")
        }

        override fun retainAll(elements: Collection<MutableMap.MutableEntry<Int, Double>>): Boolean {
            TODO("Not yet implemented")
        }

        override fun contains(element: MutableMap.MutableEntry<Int, Double>): Boolean {
            TODO("Not yet implemented")
        }

        override fun containsAll(elements: Collection<MutableMap.MutableEntry<Int, Double>>): Boolean {
            TODO("Not yet implemented")
        }

        override fun isEmpty(): Boolean {
            TODO("Not yet implemented")
        }

    }

    inline class IteratorWrapper(val apacheIterator: OpenIntToDoubleHashMap.Iterator): MutableIterator<MutableMap.MutableEntry<Int, Double>> {
        override fun hasNext() = apacheIterator.hasNext()

        override fun next(): MutableMap.MutableEntry<Int, Double> {
            return Entry(apacheIterator.key(), apacheIterator.value())
        }

        override fun remove() {
            TODO("Not yet implemented")
        }

    }

    class Entry(override val key: Int, override val value: Double): MutableMap.MutableEntry<Int,Double> {
        override fun setValue(newValue: Double): Double {
            TODO("Not yet implemented")
        }
    }
}