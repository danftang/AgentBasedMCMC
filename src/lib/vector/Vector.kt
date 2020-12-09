package lib.vector

import lib.abstractAlgebra.*
import org.apache.commons.math3.FieldElement

// Dense vector
interface Vector<T>: AlgebraicOperators<T> {
    val size: Int
    operator fun get(index: Int): T

    fun new(size: Int, init: (Int)->T): MutableVector<T> // = ArrayVector(size, init)

    operator fun plus(other: Vector<T>): Vector<T> = new(size) { i -> get(i) + other[i] }
    operator fun minus(other: Vector<T>): Vector<T> = new(size) { i -> get(i) - other[i] }
    operator fun times(multiplier: T): Vector<T> = new(size) { i -> get(i) * multiplier }

    fun asList(): List<T> {
        return object: AbstractList<T>() {
            override val size: Int
                get() = this@Vector.size
            override fun get(index: Int): T = this@Vector.get(index)
        }
    }
}


fun<T: FieldElement<T>> vectorOf(vararg values: T): Vector<T> {
    return ArrayVector(values)
}

//fun fieldVectorOf(vararg values: Double): Vector<DoubleFieldElement> {
//    return ListVector(object: AbstractList<DoubleFieldElement>() {
//        override val size: Int
//            get() = values.size
//        override fun get(index: Int) = DoubleFieldElement(values[index])
//    })
//}

fun vectorOf(vararg values: Double): Vector<Double> {
    return object: Vector<Double>, DoubleOperators {
        override val size: Int
            get() = values.size
        override fun get(index: Int) = values[index]
        override fun new(size: Int, init: (Int) -> Double) = DoubleArrayVector(size, init)
    }
}

//fun<T: Any> vectorOf(operators: FieldOperators<T>, vararg values: T): Vector<out T> {
//    return ArrayVector<T>(operators, values)
//}

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
