package lib.sparseIntMatrix

class IntArrayVector(val data: IntArray): AbstractList<Int>() {

    override val size: Int
        get() = data.size

    constructor(initialVals: List<Int>): this(IntArray(initialVals.size) {
        initialVals[it]
    })

    constructor(size: Int, initFunc: (Int) -> Int): this(IntArray(size, initFunc))

    constructor(size: Int): this(IntArray(size) { 0 })

    operator fun plus(other: IntArrayVector): IntArrayVector {
        assert(this.size == other.size)
        return IntArrayVector(this.size) {
            this[it] + other[it]
        }
    }

    operator fun plus(other: SparseIntVector): IntArrayVector {
        val result = IntArrayVector(data.copyOf())
        for(entry in other) {
            result[entry.key] += entry.value
        }
        return result
    }

    operator fun minus(other: IntArrayVector): IntArrayVector {
        assert(this.size == other.size)
        return IntArrayVector(this.size) {
            this[it] - other[it]
        }
    }

    operator fun minus(other: SparseIntVector): IntArrayVector {
        val result = IntArrayVector(data.copyOf())
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