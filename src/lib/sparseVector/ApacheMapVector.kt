package lib.sparseVector

import lib.abstractAlgebra.DoubleOperators
import org.apache.commons.math3.linear.OpenMapRealVector
import org.apache.commons.math3.util.OpenIntToDoubleHashMap

class ApacheMapVector(val apacheMap: OpenIntToDoubleHashMap): MutableSparseVector<Double>, DoubleOperators {
    override val nonZeroEntries: Map<Int,Double> = ApacheMapWrapper(apacheMap)

    constructor(): this(OpenIntToDoubleHashMap(0.0))

    constructor(vararg init: Pair<Int,Double>): this() {
        init.forEach {
            set(it.first, it.second)
        }
    }

    override fun new(): MutableSparseVector<Double>  { return ApacheMapVector(OpenIntToDoubleHashMap(0.0)) }

    override operator fun set(index: Int, value: Double) {
        if(value.isZero()) apacheMap.remove(index) else apacheMap.put(index, value)
    }

    override fun setToZero() {
        val iter = apacheMap.iterator()
        val keys = IntArray(apacheMap.size()) {
            iter.advance()
            iter.key()
        }
        keys.forEach { apacheMap.remove(it) }
    }


    override fun mapAssign(elementTransform: (Double) -> Double) {
        val iter = apacheMap.iterator()
        while(iter.hasNext()) {
            iter.advance()
            apacheMap.put(iter.key(), elementTransform(iter.value()))
        }
    }


    override fun mapAssignWithIndex(transform: (Int, Double) -> Double) {
        val iter = apacheMap.iterator()
        while(iter.hasNext()) {
            iter.advance()
            apacheMap.put(iter.key(), transform(iter.key(),iter.value()))
        }
    }


    override fun toString(): String {
        return nonZeroEntries.toString()
    }


    class ApacheMapWrapper(val apacheMap: OpenIntToDoubleHashMap): AbstractMap<Int,Double>() {
        override val entries: Set<Map.Entry<Int, Double>>
            get() = ApacheEntrySet(apacheMap)
    }


    class ApacheEntrySet(val apacheMap: OpenIntToDoubleHashMap): AbstractSet<Map.Entry<Int,Double>>() {
        override val size: Int
            get() = apacheMap.size()

        override fun iterator(): Iterator<Map.Entry<Int, Double>> {
            return IteratorWrapper(apacheMap.iterator())
        }
    }


    inline class IteratorWrapper(val apacheIterator: OpenIntToDoubleHashMap.Iterator): Iterator<Map.Entry<Int, Double>> {
        override fun hasNext() = apacheIterator.hasNext()

        override fun next(): Map.Entry<Int, Double> {
            apacheIterator.advance()
            return Entry(apacheIterator.key(), apacheIterator.value())
        }
    }


    class Entry(override val key: Int, override val value: Double): Map.Entry<Int,Double>
}