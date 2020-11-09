package lib.sparseMatrix

class IntVector(val data: IntArray): AbstractList<Int>() {

    override val size: Int
        get() = data.size

    constructor(initialVals: List<Int>): this(IntArray(initialVals.size) {
        initialVals[it]
    })

    constructor(size: Int, initFunc: (Int) -> Int): this(IntArray(size, initFunc))

    constructor(size: Int): this(IntArray(size) { 0 })

    operator fun plus(other: IntVector): IntVector {
        assert(this.size == other.size)
        return IntVector(this.size) {
            this[it] + other[it]
        }
    }

    operator fun plus(other: SparseIntVector): IntVector {
        val result = IntVector(data.copyOf())
        for(entry in other) {
            result[entry.key] += entry.value
        }
        return result
    }

    operator fun minus(other: IntVector): IntVector {
        assert(this.size == other.size)
        return IntVector(this.size) {
            this[it] - other[it]
        }
    }

    operator fun minus(other: SparseIntVector): IntVector {
        val result = IntVector(data.copyOf())
        for(entry in other) {
            result[entry.key] -= entry.value
        }
        return result
    }


    override operator fun get(index: Int): Int =
        data[index]

    operator fun set(index: Int, value: Int) {
        data[index] = value
    }

    fun toSparseIntVector(): HashIntVector {
        val result = HashIntVector()
        for(i in indices) {
            val v = this[i]
            if(v != 0) result[i] = v
        }
        return result
    }

    override fun toString(): String {
        return data.asList().toString()
    }

}