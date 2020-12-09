package lib.vector

import org.apache.commons.math3.FieldElement

class MutableArrayVector<T: FieldElement<T>>(override val data: Array<T>): ArrayVector<T>(data), MutableVector<T> {

    constructor(size: Int, init: (Int) -> T): this(Array<Any>(size, init) as Array<T>)

    override fun set(index: Int, value: T) {
        data[index] = value
    }
}