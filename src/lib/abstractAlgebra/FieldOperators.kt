package lib.abstractAlgebra

import lib.vector.MutableMapVector
import kotlin.reflect.KClass

interface FieldOperators<T>: AlgebraicOperators<T> {
    val zero: T
    val one: T
    val operators: FieldOperators<T>
        get() = this

    fun T.isZero(): Boolean  // ...in order to deal with the problem that -0.0 != 0.0 for boxed Doubles

//    val runtimeKClass: KClass<out T>

//    operator fun T.plus(other: T): T
//    operator fun T.minus(other: T): T
//    operator fun T.unaryMinus(): T
//    operator fun T.times(other: T): T
//    operator fun T.div(other: T): T
}