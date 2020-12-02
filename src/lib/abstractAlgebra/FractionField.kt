package lib.abstractAlgebra

import lib.unaryMinus
import org.apache.commons.math3.fraction.Fraction

inline class FractionField(val x: Fraction): FieldElement<FractionField> {

    override fun plus(other: FractionField)  = FractionField(x.add(other.x))
    override fun minus(other: FractionField) = FractionField(x.subtract(other.x))
    override fun unaryMinus()                = FractionField(x.unaryMinus())
    override fun times(n: Int)               = FractionField(x.multiply(n))
    override fun times(other: FractionField) = FractionField(x.multiply(other.x))
    override fun div(other: FractionField)   = FractionField(x.divide(other.x))
    override fun getField()                  = FractionField

    override fun compareTo(other: FractionField) = x.compareTo(other.x)

    fun asFraction(): Fraction = x

    companion object: lib.abstractAlgebra.Field<FractionField> {
        override fun getZero()          = FractionField(Fraction.ZERO)
        override fun getOne()           = FractionField(Fraction.ONE)
        override fun getRuntimeClass()  = FractionField::class.java
        override fun getRuntimeKClass() = FractionField::class
    }



}