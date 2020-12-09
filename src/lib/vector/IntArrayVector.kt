package lib.vector

import lib.abstractAlgebra.IntFieldElement
import lib.abstractAlgebra.IntOperators
import kotlin.reflect.KClass

inline class IntArrayVector(val data: IntArray): MutableVector<Int>, IntOperators {
    override fun set(index: Int, value: Int) {
        data[index] = value
    }

    override val size: Int
        get() = data.size

    override fun get(index: Int): Int {
        return data[index]
    }

    override fun new(size: Int, init: (Int)->Int): MutableVector<Int> {
        return IntArrayVector(IntArray(size) { init(it) })
    }
}