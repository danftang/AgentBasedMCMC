package lib.vector

import lib.abstractAlgebra.FieldOperators
import lib.abstractAlgebra.WithFieldElementOperators
import org.apache.commons.math3.FieldElement

open class ArrayVector<T: FieldElement<T>>(open val data: Array<out T>): Vector<T>, WithFieldElementOperators<T> {


    override val size: Int = data.size

    override fun get(index: Int) = data[index]
    override fun new(size: Int, init: (Int) -> T): MutableVector<T> {
        return MutableArrayVector(size, init)
    }


//    override fun new(): MutableVector<T> = ArrayVector(size)
}

inline fun<T: FieldElement<T>> Array<T>.asArrayVector(): ArrayVector<T> {
    return ArrayVector( this)
}


fun<T: FieldElement<T>> arrayVectorOf(vararg values: T): ArrayVector<T> {
    return ArrayVector(values)
}

