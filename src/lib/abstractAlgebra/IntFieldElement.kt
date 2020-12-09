package lib.abstractAlgebra

import org.apache.commons.math3.Field
import org.apache.commons.math3.FieldElement

class IntFieldElement(val value: Int): FieldElement<IntFieldElement>, Comparable<IntFieldElement>, Number() {
    override fun add(p0: IntFieldElement) = IntFieldElement(value + p0.value)
    override fun subtract(p0: IntFieldElement) = IntFieldElement(value - p0.value)
    override fun negate() = IntFieldElement(-value)
    override fun multiply(p0: Int) = IntFieldElement(value * p0)
    override fun multiply(p0: IntFieldElement) = IntFieldElement(value * p0.value)
    override fun divide(p0: IntFieldElement) = IntFieldElement(value / p0.value)
    override fun reciprocal() = IntFieldElement(1/value)

    fun asInt(): Int = value

    override fun getField(): Field<IntFieldElement> {
        TODO("Not yet implemented")
    }

    override fun compareTo(other: IntFieldElement) = value.compareTo(other.value)
    override fun equals(other: Any?): Boolean {
        return if(other is IntFieldElement) {
            value == other.value
        } else false
    }

    override fun toByte() = value.toByte()
    override fun toChar() = value.toChar()
    override fun toDouble() = value.toDouble()
    override fun toFloat() = value.toFloat()
    override fun toInt() = value
    override fun toLong() = value.toLong()
    override fun toShort() = value.toShort()

}

fun Int.asIntFieldElement(): IntFieldElement = IntFieldElement(this)
