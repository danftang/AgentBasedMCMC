package lib.abstractAlgebra

import kotlin.reflect.KClass

object DoubleOperators: FieldOperators<Double> {
    override fun Double.plus(other: Double) = this + other
    override fun Double.minus(other: Double) = this - other
    override fun Double.unaryMinus() = -this
    override fun Double.times(other: Double) = this * other
    override fun Double.div(other: Double) = this / other
    override fun Double.compareTo(other: Double) = this.compareTo(other)

    override val zero: Double
        get() = 0.0

    override val one: Double
        get() = 1.0

    override val runtimeClass: KClass<Double>
        get() = Double::class
}