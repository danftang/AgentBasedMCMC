package lib

import java.lang.RuntimeException

inline fun Number.plus(other: Number) : Number = when(this) {
        is Double   -> this + other.toDouble()
        is Int      -> this + other.toInt()
        else        -> throw(RuntimeException("Don't know how to add these things"))
    }

interface ArithmeticContext<T> {
    infix fun T.plus(b: T): T
}

fun<T> test(context: ArithmeticContext<T>, a: T, b: T) {
    with(context) {
        val c = a plus b
    }
}