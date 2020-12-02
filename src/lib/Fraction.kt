package lib

import org.apache.commons.math3.fraction.Fraction

// Extensions for Apache commons math Fraction class

operator fun Fraction.plus(i: Int) = this.add(i)

operator fun Fraction.minus(i: Int) = this.subtract(i)

operator fun Fraction.plus(i: Fraction) = this.add(i)

operator fun Fraction.minus(i: Fraction) = this.subtract(i)

operator fun Fraction.times(i: Int) = this.multiply(i)

operator fun Fraction.times(other: Fraction) = this.multiply(other)

operator fun Fraction.div(i: Int) = this.divide(i)

operator fun Fraction.unaryMinus() = this.negate()

fun Fraction.isInteger(): Boolean = numerator.rem(denominator) == 0

fun max(a: Fraction, b: Fraction) = if(b > a) b else a
