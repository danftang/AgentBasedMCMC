package lib.abstractAlgebra

import kotlin.reflect.KClass

interface IntOperators: FieldOperators<Int> {
    override fun Int.plus(other: Int) = this + other
    override fun Int.minus(other: Int) = this - other
    override fun Int.unaryMinus() = -this
    override fun Int.times(other: Int) = this * other
    override fun Int.div(other: Int) = this / other
    override fun Int.isZero() = this == 0

    override val zero: Int
        get() = 0

    override val one: Int
        get() = 1

    override val operators
        get() = IntOperators

    override fun Int.toDouble(): Double = this.toDouble()
    override fun Int.toInt(): Int = this
    override fun Int.toField(): Int = this
    override fun Double.toField(): Int = this.toInt()

    companion object: IntOperators
}