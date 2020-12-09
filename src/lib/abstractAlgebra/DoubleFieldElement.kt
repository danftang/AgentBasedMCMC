package lib.abstractAlgebra

import org.apache.commons.math3.Field
import org.apache.commons.math3.FieldElement

class DoubleFieldElement(val value: Double): FieldElement<DoubleFieldElement>, Comparable<DoubleFieldElement>, Number() {
    override fun add(p0: DoubleFieldElement) = DoubleFieldElement(value + p0.value)
    override fun subtract(p0: DoubleFieldElement) = DoubleFieldElement(value - p0.value)
    override fun negate() = DoubleFieldElement(-value)
    override fun multiply(p0: Int) = DoubleFieldElement(value * p0)
    override fun multiply(p0: DoubleFieldElement) = DoubleFieldElement(value * p0.value)
    override fun divide(p0: DoubleFieldElement) = DoubleFieldElement(value / p0.value)
    override fun reciprocal() = DoubleFieldElement(1.0/value)

    fun asDouble(): Double = value

    override fun getField(): Field<DoubleFieldElement> {
        TODO("Not yet implemented")
    }

    override fun compareTo(other: DoubleFieldElement) = value.compareTo(other.value)
    override fun equals(other: Any?): Boolean {
        return if(other is DoubleFieldElement) {
            value == other.value
        } else false
    }

    override fun toByte() = value.toByte()
    override fun toChar() = value.toChar()
    override fun toDouble() = value
    override fun toFloat() = value.toFloat()
    override fun toInt() = value.toInt()
    override fun toLong() = value.toLong()
    override fun toShort() = value.toShort()
}

fun Double.asDoubleFieldElement(): DoubleFieldElement = DoubleFieldElement(this)
