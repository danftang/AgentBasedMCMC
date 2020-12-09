package lib.abstractAlgebra

import org.apache.commons.math3.FieldElement

interface WithFieldElementOperators<T: FieldElement<T>>: AlgebraicOperators<T> {
    override fun T.plus(other: T): T = this.add(other)
    override fun T.minus(other: T): T = this.subtract(other)
    override fun T.unaryMinus(): T = this.negate()
    override fun T.times(other: T): T = this.multiply(other)
    override fun T.div(other: T): T = this.divide(other)
}