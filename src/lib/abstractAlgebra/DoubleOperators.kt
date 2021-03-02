package lib.abstractAlgebra

import kotlin.reflect.KClass

interface DoubleOperators: FieldOperators<Double> {
    override fun Double.plus(other: Double) = this + other
    override fun Double.minus(other: Double) = this - other
    override fun Double.unaryMinus() = -this
    override fun Double.times(other: Double) = this * other
    override fun Double.div(other: Double) = this / other
    override fun Double.isZero() = this == 0.0

    override val zero: Double
        get() = 0.0

    override val one: Double
        get() = 1.0

    override val operators
        get() = DoubleOperators

    override fun Double.toDouble(): Double = this
    override fun Double.toInt(): Int = this.toInt()
    override fun Int.toField(): Double = this.toDouble()
    override fun Double.toField(): Double = this

    companion object: DoubleOperators
}