package lib.vector

import lib.abstractAlgebra.DoubleOperators

inline class DoubleArrayVector(val data: DoubleArray): MutableVector<Double>, DoubleOperators {

    constructor(size: Int, init: (Int)->Double): this(DoubleArray(size, init))

    override fun set(index: Int, value: Double) {
        data[index] = value
    }

    override val size: Int
        get() = data.size

    override fun get(index: Int) = data[index]

    override fun new(size: Int, init: (Int) -> Double): MutableVector<Double> {
        return DoubleArrayVector(size, init)
    }
}