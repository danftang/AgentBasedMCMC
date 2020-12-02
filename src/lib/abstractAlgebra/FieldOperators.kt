package lib.abstractAlgebra

import kotlin.reflect.KClass

interface FieldOperators<T: Any> {
    val zero: T
    val one: T
    val runtimeClass: KClass<T>

    operator fun T.plus(other: T): T
    operator fun T.minus(other: T): T
    operator fun T.unaryMinus(): T
    operator fun T.times(other: T): T
    operator fun T.div(other: T): T
    operator fun T.compareTo(other: T): Int

}