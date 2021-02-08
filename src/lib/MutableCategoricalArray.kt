import kotlin.random.Random

// Represents a probability distribution over an integer range 0..N
// Can be treated as an array of doubles, where each double is an
// (un-normalised) probability for that index.
//
// Internally this is stored as a binary sum tree. However, we
// only store sums of nodes that are right hand children. This allows
// us to map the tree onto a DoubleArray. The bits of the index of a node
// can be calculated, starting from the most significant bit, from the sequence of
// left/right (0/1) branches on the path from the root to that node.
// Since the path must end with a right branch (1) we pad the sequence with
// zeroes beyond the final right branch to give a unique index.
//
// This encoding allows arrays of any size (not just integer multiples of 2)
// and allows modification of probabilities and sampling in O(log(N)) time.
// If all probabilities need modifying, this can be done in O(N) time using
// the setAll method.
class MutableCategoricalArray: AbstractList<Double> {

    val random: Random
    override val size: Int
        get() = tree.size

    private val tree: DoubleArray
    private val indexHighestBit: Int


    constructor(size: Int, random: Random = Random) {
        tree = DoubleArray(size)
        this.random = random
        indexHighestBit = (tree.size-1).highestOneBit()
    }


    constructor(size: Int, init: (Int) -> Double): this(size, Random, init)


    constructor(size: Int, random: Random, init: (Int) -> Double): this(size,random) {
        for(i in size-1 downTo 0) { tree[i] = descendantSum(i) + init(i) }
    }


    // sets the un-normalised probability associated with the supplied index
    operator fun set(index: Int, probability: Double) {
        var sum = probability
        var indexOffset = 1
        while((indexOffset and index) == 0 && indexOffset < size) {
            val descendantIndex = index + indexOffset
            if(descendantIndex < size) sum += tree[descendantIndex]
            indexOffset = indexOffset shl 1
        }
        val delta = sum - tree[index]
        var ancestorIndex = index
        tree[index] = sum
        while(indexOffset < size) {
            ancestorIndex = ancestorIndex xor indexOffset
            tree[ancestorIndex] += delta
            do {
                indexOffset = indexOffset shl 1
            } while(ancestorIndex and indexOffset == 0 && indexOffset < size)
        }
    }


    // returns the un-normalised probability associated with the supplied index.
    override operator fun get(index: Int): Double = tree[index] - descendantSum(index)


    // draws a sample from the distribution
    fun sample(): Int {
        var index = 0
        var target = random.nextDouble() * tree[0]
        var rightChildOffset = indexHighestBit
        while(rightChildOffset != 0) {
            val childIndex = index+rightChildOffset
            if(childIndex < size) {
                if (tree[childIndex] > target) index += rightChildOffset else target -= tree[childIndex]
            }
            rightChildOffset = rightChildOffset shr 1
        }
        return index
    }


    // Sets the un-normalised probabilities of the first N integers
    fun setAll(values: List<Double>) {
        for(i in values.size-1 downTo 0) {
            tree[i] = descendantSum(i) + values[i]
        }
    }


    // the sum of all un-normalised probabilities (doesn't need to be 1.0)
    fun sum(): Double {
        return tree[0]
    }


    // Returns the normalised probability of the index'th element
    fun P(index: Int): Double = get(index)/sum()


//    fun asList(): List<Double> {
//        return object: AbstractList<Double>() {
//            override val size: Int
//                get() = this@MutableCategoricalArray.size
//            override fun get(index: Int) = this@MutableCategoricalArray[index]
//        }
//    }


    private fun descendantSum(index: Int): Double {
        var indexOffset = 1
        var sum = 0.0
        while((indexOffset and index) == 0 && indexOffset < size) {
            val descendantIndex = index + indexOffset
            if(descendantIndex < size) sum += tree[descendantIndex]
            indexOffset = indexOffset shl 1
        }
        return sum
    }


    private fun Int.highestOneBit(): Int {
        var i = this
        i = i or (i shr 1)
        i = i or (i shr 2)
        i = i or (i shr 4)
        i = i or (i shr 8)
        i = i or (i shr 16)
        return i - (i ushr 1)
    }
}

fun mutableCategoricalOf(vararg probabilities: Double): MutableCategoricalArray {
    return MutableCategoricalArray(probabilities.size) { probabilities[it] }
}
