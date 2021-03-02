package lib.abstractAlgebra

import org.apache.commons.math3.fraction.Fraction

interface FractionOperators: FieldOperators<Fraction> {

    override val zero: Fraction
        get() = Fraction.ZERO
        override val one: Fraction
          get() = Fraction.ONE
       override val operators: FieldOperators<Fraction>
         get() = this

    override fun Fraction.plus(other: Fraction) = this.add(other)
    override fun Fraction.minus(other: Fraction) = this.subtract(other)
    override fun Fraction.unaryMinus() = this.negate()
    override fun Fraction.times(other: Fraction) = this.multiply(other)
    override fun Fraction.div(other: Fraction) = this.divide(other)
    override fun Fraction.isZero() = this == zero

    override fun Fraction.toDouble(): Double = this.toDouble()
    override fun Fraction.toInt(): Int = this.toInt()
    override fun Int.toField(): Fraction = Fraction(this)
    override fun Double.toField(): Fraction = Fraction(this)

    companion object: FractionOperators
}