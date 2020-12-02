package lib.abstractAlgebra

import kotlin.reflect.KClass

object IntOperators: FieldOperators<Int> {
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