package lib.vector

import lib.abstractAlgebra.Field
import lib.abstractAlgebra.FieldElement
import lib.abstractAlgebra.FieldOperators

class ArrayVector<T: Any>(val fieldOperators: FieldOperators<T>, val data: Array<T>): MutableVector<T>, FieldOperators<T> by fieldOperators {

    constructor(fieldOperators: FieldOperators<T>, size: Int): this(fieldOperators, Array<Any>(size) { fieldOperators.zero } as Array<T>)

    override val size: Int = data.size

    override fun get(index: Int) = data[index]

    override fun set(index: Int, value: T) {
        data[index] = value
    }

    override fun new(): MutableVector<T> = ArrayVector(fieldOperators, size)
}