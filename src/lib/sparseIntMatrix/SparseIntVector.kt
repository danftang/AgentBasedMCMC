package lib.sparseIntMatrix

interface SparseIntVector: Iterable<Map.Entry<Int,Int>> {
    val sparseSize: Int     // number of non-zero entries

    val keys: Iterable<Int>
        get() = this.asSequence().map { it.key }.asIterable()

    val values: Iterable<Int>
        get() = this.asSequence().map { it.value }.asIterable()


    operator fun get(index: Int): Int

    fun toIntArray(size: Int): IntArray {
        val array = IntArray(size)
        for(entry in this) {
            array[entry.key] = entry.value
        }
        return array
    }


    fun mapKeys(mapFunction: (Map.Entry<Int, Int>) -> Int): HashIntVector {
        val mappedVector = HashIntVector(sparseSize)
        for(entry in this) {
            mappedVector[mapFunction(entry)] = entry.value
        }
        return mappedVector
    }

    fun mapValues(mapFunction: (Map.Entry<Int,Int>) -> Int): HashIntVector {
        val mappedVector = HashIntVector(sparseSize)
        for(entry in this) {
            mappedVector[entry.key] = mapFunction(entry)
        }
        return mappedVector
    }

    operator fun plus(other: SparseIntVector): HashIntVector {
        val sum = HashIntVector(this)
        for(entry in other) {
            sum.data.merge(entry.key, entry.value) {a,b ->
                val result = a+b
                if(result !=0) result else null
            }
        }
        return sum
    }

    operator fun plus(other: IntArrayVector): IntArrayVector {
        return other + this
    }


    operator fun minus(other: SparseIntVector): HashIntVector {
        val sum = HashIntVector(this)
        sum -= other
        return sum
    }


    operator fun minus(other: IntArrayVector): IntArrayVector {
        val result = IntArrayVector(other.size) {
            -other[it]
        }
        for(entry in this) {
            result[entry.key] += entry.value
        }
        return result
    }


    operator fun unaryMinus(): HashIntVector {
        val result = HashIntVector(HashMap(this.sparseSize))
        for (entry in this) {
            result[entry.key] = -entry.value
        }
        return result
    }



    fun isPositive(): Boolean {
        return values.fold(true) { isPositive, v -> isPositive && (v >= 0) }
    }

    fun normSqL2(): Int = this.values.sumBy { it*it }

    fun dotProd(other: SparseIntVector): Int {
        val columns =
            if(this.sparseSize < other.sparseSize) Pair(this, other) else Pair(other,this)
        return columns.first.sumBy { it.value * columns.second[it.key] }
    }

    fun dotProd(other: IntArrayVector): Int {
        var sum = 0
        for(i in other.indices) {
            sum += other[i]*this[i]
        }
        return sum
    }


    fun valueRange(): IntRange {
        var min = Int.MAX_VALUE
        var max = Int.MIN_VALUE
        values.forEach {
            if(it > max) max = it
            if(it < min) min = it
        }
        return IntRange(min, max)
    }

    fun toIntVector(size: Int): IntArrayVector {
        val result = IntArrayVector(size)
        for(entry in this) {
            result[entry.key] = entry.value
        }
        return result
    }

}