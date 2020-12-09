package lib.abstractAlgebra

interface AlgebraicOperators<T> {
    operator fun T.plus(other: T): T
    operator fun T.minus(other: T): T
    operator fun T.unaryMinus(): T
    operator fun T.times(other: T): T
    operator fun T.div(other: T): T
}