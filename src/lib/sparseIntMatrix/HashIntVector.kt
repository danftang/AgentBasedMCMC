package lib.sparseIntMatrix

class HashIntVector(val data: HashMap<Int,Int>): SparseIntVector {
    override val sparseSize: Int
        get() = data.size

    constructor() : this(HashMap(4))

    constructor(initialCapacity: Int) : this(HashMap(initialCapacity))

    constructor(copy: SparseIntVector) : this(HashMap(copy.sparseSize)) {
        for (entry in copy) {
            data[entry.key] = entry.value
        }
    }

    constructor(copy: IntArrayVector) : this() {
        for (i in copy.indices) {
            val xi = copy[i]
            if(xi != 0) data[i] = xi
        }
    }

    constructor(vararg entries: Pair<Int, Int>) : this(HashMap(entries.size)) {
        data.putAll(entries)
    }

    override fun get(index: Int): Int = data[index] ?: 0

    operator fun set(index: Int, value: Int) {
        if(value == 0) {
            data.remove(index)
            return
        }
        data[index] = value
    }

    fun setEqual(denseVector: IntArrayVector) {
        data.clear()
        for(i in 0 until denseVector.size) {
            val xi = denseVector[i]
            if(xi != 0) data[i] = xi
        }
    }

    override fun iterator(): Iterator<Map.Entry<Int, Int>> = data.entries.iterator()

    fun weightedPlusAssign(other: SparseIntVector, weight: Int) {
        for (entry in other) {
            data.merge(entry.key, entry.value * weight) { a, b ->
                val result = a + b
                if (result != 0) result else null
            }
        }
    }

    fun clear() {
        data.clear()
    }

    override fun equals(other: Any?): Boolean {
        return if (other is SparseIntVector) {
            this.sparseSize == other.sparseSize && this.all { other[it.key] == it.value }
        } else false
    }

    operator fun timesAssign(multiplier: Int) {
        if (multiplier == 0) data.clear()
        if (multiplier == 1) return
        for (entry in data) {
            entry.setValue(entry.value * multiplier)
        }
    }

    operator fun divAssign(denominator: Int) {
        if (denominator == 1) return
        for (entry in data) {
            entry.setValue(entry.value / denominator)
        }
    }

    operator fun minusAssign(otherCol: SparseIntVector) {
        for (entry in otherCol) {
            data.merge(entry.key, -entry.value) { a, b ->
                val result = a + b
                if (result != 0) result else null
            }
        }
    }

    operator fun plusAssign(otherCol: SparseIntVector) {
        for (entry in otherCol) {
            data.merge(entry.key, entry.value) { a, b ->
                val result = a + b
                if (result != 0) result else null
            }
        }
    }

    override fun toString(): String {
        return data.toString()
    }

}