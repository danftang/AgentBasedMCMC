package lib.abstractAlgebra

interface FieldOperators<T>: AlgebraicOperators<T> {
    val zero: T
    val one: T
    val operators: FieldOperators<T>
        get() = this

    fun T.isZero(): Boolean  // ...in order to deal with the problem that -0.0 != 0.0 for boxed Doubles

    fun T.toDouble(): Double
    fun T.toInt(): Int
    fun Int.toField(): T
    fun Double.toField(): T

//    val runtimeKClass: KClass<out T>

}