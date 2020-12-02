package lib.vector

import lib.abstractAlgebra.Field
import lib.abstractAlgebra.FieldElement
import lib.abstractAlgebra.FieldOperators

// Dense vector
interface Vector<T: Any>: FieldOperators<T> {
    val size: Int
    operator fun get(index: Int): T

    fun new(): MutableVector<T>

    operator fun plus(other: Vector<T>): Vector<T> {
        val result = new()
        for (i in 0 until size) result[i] = get(i) + other[i]
        return result
    }

    operator fun minus(other: Vector<T>): Vector<T> {
        val result = new()
        for(i in 0 until size) result[i] = this[i] - other[i]
        return result
    }

    operator fun times(multiplier: T): Vector<T> {
        val result = new()
        for(i in 0 until size) result[i] = this[i] * multiplier
        return result
    }

    fun asList(): List<T> {
        return object: AbstractList<T>() {
            override val size: Int
                get() = this@Vector.size
            override fun get(index: Int): T = this@Vector.get(index)
        }
    }
}


//inline fun<FIELD, T> Vector<FIELD,T>.map(other: Vector<FIELD,T>, transform: (T,T)->T): MutableVector<FIELD,T> {
//    val result = new()
//    for(i in 0 until size) {
//        result[i] = transform(this[i], other[i])
//    }
//    return result
//}
//
//inline fun<T> Vector<FIELD,T>.map(scalar: T, transform: (T,T)->T): MutableVector<FIELD,T> {
//    val result = new()
//    for(i in 0 until size) {
//        result[i] = transform(this[i], scalar)
//    }
//    return result
//}
