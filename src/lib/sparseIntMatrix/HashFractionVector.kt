package lib.sparseIntMatrix

import lib.isInteger
import lib.sparseIntMatrix.HashIntVector
import lib.sparseIntMatrix.SparseIntVector
import org.apache.commons.math3.fraction.Fraction
import kotlin.math.roundToInt

class HashFractionVector(val data: HashMap<Int, Fraction>): Iterable<Map.Entry<Int,Fraction>> {
    val sparseSize: Int
        get() = data.size

    val keys: Iterable<Int>
        get() = this.asSequence().map { it.key }.asIterable()

    val values: Iterable<Fraction>
        get() = this.asSequence().map { it.value }.asIterable()



    constructor() : this(HashMap(4))

    constructor(copy: SparseIntVector): this(
        copy.associateByTo(HashMap<Int,Fraction>(), { it.key }, { Fraction(it.value,1) })
    )

    constructor(initialCapacity: Int) : this(HashMap(initialCapacity))

    constructor(vararg entries: Pair<Int, Fraction>) : this(HashMap(entries.size)) {
        data.putAll(entries)
    }

    operator fun get(index: Int): Fraction = data[index] ?: Fraction.ZERO

    operator fun set(index: Int, value: Fraction) {
        if(value == Fraction.ZERO) {
            data.remove(index)
            return
        }
        data[index] = value
    }

    operator fun set(index: Int, value: Int) {
        if(value == 0) {
            data.remove(index)
            return
        }
        data[index] = Fraction(value,1)
    }

    override fun iterator(): Iterator<Map.Entry<Int, Fraction>> = data.entries.iterator()

//    fun weightedPlusAssign(other: SparseIntVector, weight: Int) {
//        for (entry in other) {
//            data.merge(entry.key, entry.value * weight) { a, b ->
//                val result = a + b
//                if (result != 0) result else null
//            }
//        }
//    }

    fun roundToIntVector(): HashIntVector {
        val intVec = HashIntVector(this.sparseSize)
        for((i, xi) in this) {
            intVec[i] = xi.toDouble().roundToInt()
        }
        return intVec
    }

    fun isPositive(): Boolean {
        return all { it.value >= Fraction.ZERO }
    }

    fun isInteger(): Boolean {
        return all { it.value.isInteger() }
    }

    fun clear() {
        data.clear()
    }

    override fun equals(other: Any?): Boolean {
        return if (other is HashFractionVector) {
            this.sparseSize == other.sparseSize && this.all { other[it.key] == it.value }
        } else false
    }

//    operator fun timesAssign(multiplier: Int) {
//        if (multiplier == 0) data.clear()
//        if (multiplier == 1) return
//        for (entry in data) {
//            entry.setValue(entry.value * multiplier)
//        }
//    }
//
//    operator fun divAssign(denominator: Int) {
//        if (denominator == 1) return
//        for (entry in data) {
//            entry.setValue(entry.value / denominator)
//        }
//    }
//
//    operator fun minusAssign(otherCol: SparseIntVector) {
//        for (entry in otherCol) {
//            data.merge(entry.key, -entry.value) { a, b ->
//                val result = a + b
//                if (result != 0) result else null
//            }
//        }
//    }
//
//    operator fun plusAssign(otherCol: SparseIntVector) {
//        for (entry in otherCol) {
//            data.merge(entry.key, entry.value) { a, b ->
//                val result = a + b
//                if (result != 0) result else null
//            }
//        }
//    }

    override fun toString(): String {
        return data.toString()
    }

}