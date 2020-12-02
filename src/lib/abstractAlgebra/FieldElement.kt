package lib.abstractAlgebra

import org.apache.commons.math3.Field

interface FieldElement<T: FieldElement<T>>: org.apache.commons.math3.FieldElement<T>, Comparable<T> {
    operator fun plus(other: T): T
    operator fun minus(other: T): T
    operator fun unaryMinus(): T
    operator fun times(n: Int): T
    operator fun times(other: T): T
    operator fun div(other: T): T
    override fun getField(): Field<T>

    override fun add(a: T): T = this + a
    override fun subtract(a: T): T = this - a
    override fun negate(): T = -this
    override fun multiply(n: Int): T = this * n
    override fun multiply(a: T): T = this * a
    override fun divide(a: T): T = this / a
    override fun reciprocal(): T = field.one / (this as T)

}