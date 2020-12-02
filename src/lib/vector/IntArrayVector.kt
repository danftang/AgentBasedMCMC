package lib.vector

import kotlin.reflect.KClass

inline class IntArrayVector(val data: IntArray): MutableVector<Int> {
    override fun set(index: Int, value: Int) {
        data[index] = value
    }

    override val size: Int
        get() = data.size

    override fun get(index: Int): Int {
        return data[index]
    }

    override fun new(): MutableVector<Int> {
        return IntArrayVector(IntArray(size) { 0 })
    }

    // FieldOperators<Int> implementation
    override fun Int.plus(other: Int) = this + other
    override fun Int.minus(other: Int) = this - other
    override fun Int.unaryMinus() = -this
    override fun Int.times(other: Int) = this * other
    override fun Int.div(other: Int) = this / other
    override fun Int.compareTo(other: Int) = this.compareTo(other)

    override val zero: Int
        get() = 0

    override val one: Int
        get() = 1

    override val runtimeClass: KClass<Int>
        get() = Int::class

}