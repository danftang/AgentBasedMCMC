package lib.vector

import lib.abstractAlgebra.DoubleFieldElement
import org.apache.commons.math3.FieldElement

interface MutableVector<T>: Vector<T> {
    operator fun set(index: Int, value: T)

    operator fun timesAssign(multiplier: T) {
        mapAssign { it * multiplier }
    }

    operator fun divAssign(denominator: T) {
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

inline operator fun MutableVector<DoubleFieldElement>.set(index: Int, value: Double) {
    set(index, DoubleFieldElement(value))
}

inline fun<T> MutableVector<T>.mapAssign(unaryOp: (T)->T) {
    for (i in 0 until size) set(i, unaryOp(get(i)))
}

inline fun<T> MutableVector<T>.mapAssign(other: Vector<T>, binaryOp: (T, T)->T) {
    for(i in 0 until size) set(i, binaryOp(get(i),other[i]))
}
