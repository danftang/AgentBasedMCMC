package lib

import java.util.ArrayDeque

// Collection of objects of type T. Each object is
// associated with an index that is assigned by the
// collection at creation by calling the "new" method.
// Objects can be retrieved by suppying the [] operator
// with the object's index.
//
// Use this container when you have a very large number
// of objects. This will be faster than Java's memory
// management.
class ArrayStore<T>(val data: ArrayList<T>) {
    val freeIndices = ArrayDeque<Int>()

    fun new(initVal: T): Int {
        if(freeIndices.isNotEmpty()) {
            val index = freeIndices.poll()
            data[index] = initVal
            return index
        }
        data.add(initVal)
        return data.lastIndex
    }

    fun delete(index: Int) {
        freeIndices.add(index)
    }

    operator fun get(index: Int): T = data[index]
}