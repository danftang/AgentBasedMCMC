package lib.abstractAlgebra

import org.apache.commons.math3.Field
import org.apache.commons.math3.FieldElement
import org.apache.commons.math3.fraction.Fraction
import kotlin.reflect.KClass

inline class FieldElementOperators<T: FieldElement<T>>(val apacheField: Field<T>): FieldOperators<T> {
    override val zero: T
        get() = apacheField.zero
    override val one: T
        get() = apacheField.one
    override val operators: FieldOperators<T>
        get() = FieldElementOperators(apacheField)

    override fun T.plus(other: T) = this.add(other)
    override fun T.minus(other: T) = this.subtract(other)
    override fun T.unaryMinus() = this.negate()
    override fun T.times(other: T) = this.multiply(other)
    override fun T.div(other: T) = this.divide(other)
    override fun T.isZero() = this == zero
}

fun<T: FieldElement<T>> Field<T>.asFieldOperators() = FieldElementOperators(this)

val FractionOperators = Fraction.ZERO.field.asFieldOperators()
