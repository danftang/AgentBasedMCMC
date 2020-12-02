package lib.vector

import lib.abstractAlgebra.Field
import lib.abstractAlgebra.FieldElement
import lib.abstractAlgebra.FieldOperators

interface MutableVector<T: Any>: Vector<T> {
    operator fun set(index: Int, value: T)

    operator fun timesAssign(multiplier: T) {
        if (multiplier == one) return
        if (multiplier == zero) {
            for(i in 0 until size) set(i,zero)
            return
        }
        mapAssign { it * multiplier }
    }

    operator fun divAssign(denominator: T) {
        if(denominator == one) return
        mapAssign { it / denominator }
    }

    operator fun plusAssign(other: Vector<T>) {
        mapAssign(other) {a,b -> a + b}
    }

    operator fun minusAssign(other: Vector<T>) {
        mapAssign(other) {a,b -> a - b}
    }

    fun weightedPlusAssign(other: Vector<T>, weight: T) {
        mapAssign(other) {a,b -> a + weight*b}
    }

}

inline fun<T: Any> MutableVector<T>.mapAssign(unaryOp: (T)->T) {
    for (i in 0 until size) set(i, unaryOp(get(i)))
}

inline fun<T: Any> MutableVector<T>.mapAssign(other: Vector<T>, binaryOp: (T, T)->T) {
    for(i in 0 until size) set(i, binaryOp(get(i),other[i]))
}
