package lib

import java.util.*

/*
 ** KDTree.java by Julian Kent
 ** Converted to Kotlin by Daniel Tang, 2020
 **
 ** Licenced under the  Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 **
 ** Licence summary:
 ** Under this licence you are free to:
 **      Share : copy and redistribute the material in any medium or format
 **      Adapt : remix, transform, and build upon the material
 **      The licensor cannot revoke these freedoms as long as you follow the license terms.
 **
 ** Under the following terms:
 **      Attribution:
 **            You must give appropriate credit, provide a link to the license, and indicate
 **            if changes were made. You may do so in any reasonable manner, but not in any
 **            way that suggests the licensor endorses you or your use.
 **      NonCommercial:
 **            You may not use the material for commercial purposes.
 **      ShareAlike:
 **            If you remix, transform, or build upon the material, you must distribute your
 **            contributions under the same license as the original.
 **      No additional restrictions:
 **            You may not apply legal terms or technological measures that legally restrict
 **            others from doing anything the license permits.
 **
 ** See full licencing details here: http://creativecommons.org/licenses/by-nc-sa/3.0/
 **
 ** For additional licencing rights (including commercial) please contact jkflying@gmail.com
 **
 */
abstract class KDTree<T> private constructor(private val _dimensions: Int) {
    private var _nodes = 0
    private val root: Node
    private val nodeList = ArrayList<Node>()
    // prevent GC from having to collect _bucketSize*dimensions*sizeof(double) bytes each
// time a leaf splits
    private var mem_recycle: DoubleArray
    // the starting values for bounding boxes, for easy access
    private val bounds_template: DoubleArray
    // one big self-expanding array to keep all the node bounding boxes so that
// they stay in cache
// node bounds available at:
// low: 2 * _dimensions * node.index + 2 * dim
// high: 2 * _dimensions * node.index + 2 * dim + 1
    private val nodeMinMaxBounds: ContiguousDoubleArrayList

    val values: Sequence<T>
        get() = nodeList.asSequence().flatMap { it.pointPayloads?.asSequence()?:emptySequence() }

    val keys: Sequence<List<Double>>
        get() = nodeList
            .asSequence()
            .flatMap { it.pointLocations
                ?.array
                ?.asSequence()
                ?.take(size *_dimensions)
                ?.windowed(_dimensions, _dimensions)
                ?: emptySequence()
            }

    val entries: Sequence<Pair<List<Double>,T>>
        get() = keys.zip(values)

    fun nodes(): Int {
        return _nodes
    }

    val size: Int
        get() {
            return root.entries
        }

    operator fun set(vararg location: Double, value: T) {
        addPoint(location, value)
    }

    fun addPoint(location: DoubleArray, payload: T): Int {
        var addNode = root
        // Do a Depth First Search to find the Node where 'location' should be
// stored
        while (addNode.pointLocations == null) {
            addNode.expandBounds(location)
            addNode =
                if (location[addNode.splitDim] < addNode.splitVal) nodeList[addNode.lessIndex] else nodeList[addNode.moreIndex]
        }
        addNode.expandBounds(location)
        val nodeSize = addNode.add(location, payload)
        if (nodeSize % _bucketSize == 0) // try splitting again once every time the node passes a _bucketSize
// multiple
// in case it is full of points of the same location and won't split
            addNode.split()
        return root.entries
    }

    fun nearestNeighbours(K: Int, vararg searchLocation: Double): List<SearchResult<T?>> {
        var K = K
        K = Math.min(K, size)
        val returnResults =
            ArrayList<SearchResult<T?>>(K)
        if (K > 0) {
            val stack = IntStack()
            val results = PrioQueue<T>(K)
            stack.push(root.index)
            var added = 0
            while (stack.size() > 0) {
                val nodeIndex = stack.pop()
                if (added < K || results.peekPrio() > pointRectDist(nodeIndex, searchLocation)) {
                    val node = nodeList[nodeIndex]
                    if (node.pointLocations == null) node.search(searchLocation, stack) else added += node.search(
                        searchLocation,
                        results
                    )
                }
            }
            val priorities = results.priorities
            val elements = results.elements
            for (i in 0 until K) { // forward (closest first)
                val s =
                    SearchResult(priorities[i], elements[i] as T?)
                returnResults.add(s)
            }
        }
        return returnResults
    }

    fun nearestNeighbour(vararg searchLocation: Double): SearchResult<T?> {
        return nearestNeighbours(1, *searchLocation).first()
    }

    fun ballSearch(searchLocation: DoubleArray, radius: Double): ArrayList<T> {
        val stack = IntStack()
        val results = ArrayList<T>()
        stack.push(root.index)
        while (stack.size() > 0) {
            val nodeIndex = stack.pop()
            if (radius > pointRectDist(nodeIndex, searchLocation)) {
                val node = nodeList[nodeIndex]
                if (node.pointLocations == null) stack.push(node.moreIndex).push(node.lessIndex) else node.searchBall(
                    searchLocation,
                    radius,
                    results
                )
            }
        }
        return results
    }

    fun rectSearch(mins: DoubleArray, maxs: DoubleArray): ArrayList<T> {
        val stack = IntStack()
        val results = ArrayList<T>()
        stack.push(root.index)
        while (stack.size() > 0) {
            val nodeIndex = stack.pop()
            if (overlaps(mins, maxs, nodeIndex)) {
                val node = nodeList[nodeIndex]
                if (node.pointLocations == null) stack.push(node.moreIndex).push(node.lessIndex) else node.searchRect(
                    mins,
                    maxs,
                    results
                )
            }
        }
        return results
    }

    abstract fun pointRectDist(offset: Int, location: DoubleArray): Double
    abstract fun pointDist(arr: DoubleArray, location: DoubleArray, index: Int): Double

    fun contains(arr: DoubleArray, mins: DoubleArray, maxs: DoubleArray, index: Int): Boolean {
        var offset = (index + 1) * mins.size
        var i = mins.size
        while (i-- > 0) {
            val d = arr[--offset]
            if (mins[i] > d || d > maxs[i]) return false
        }
        return true
    }

    fun overlaps(mins: DoubleArray, maxs: DoubleArray, offset: Int): Boolean {
        var offset = offset
        offset *= 2 * maxs.size
        val array = nodeMinMaxBounds.array
        var i = 0
        while (i < maxs.size) {
            val bmin = array[offset]
            val bmax = array[offset + 1]
            if (mins[i] > bmax || maxs[i] < bmin) return false
            i++
            offset += 2
        }
        return true
    }

    class Euclidean<T>(dims: Int) : KDTree<T>(dims) {
        override fun pointRectDist(offset: Int, location: DoubleArray): Double {
            var offset = offset
            offset *= 2 * super._dimensions
            var distance = 0.0
            val array = super.nodeMinMaxBounds.array
            var i = 0
            while (i < location.size) {
                var diff = 0.0
                var bv = array[offset]
                val lv = location[i]
                if (bv > lv) diff = bv - lv else {
                    bv = array[offset + 1]
                    if (lv > bv) diff = lv - bv
                }
                distance += sqr(diff)
                i++
                offset += 2
            }
            return distance
        }

        override fun pointDist(arr: DoubleArray, location: DoubleArray, index: Int): Double {
            var distance = 0.0
            var offset = (index + 1) * super._dimensions
            var i = super._dimensions
            while (i-- > 0) {
                distance += sqr(arr[--offset] - location[i])
            }
            return distance
        }
    }

    class Manhattan<T>(dims: Int) : KDTree<T>(dims) {
        override fun pointRectDist(offset: Int, location: DoubleArray): Double {
            var offset = offset
            offset *= 2 * super._dimensions
            var distance = 0.0
            val array = super.nodeMinMaxBounds.array
            var i = 0
            while (i < location.size) {
                var diff = 0.0
                var bv = array[offset]
                val lv = location[i]
                if (bv > lv) diff = bv - lv else {
                    bv = array[offset + 1]
                    if (lv > bv) diff = lv - bv
                }
                distance += diff
                i++
                offset += 2
            }
            return distance
        }

        override fun pointDist(arr: DoubleArray, location: DoubleArray, index: Int): Double {
            var distance = 0.0
            var offset = (index + 1) * super._dimensions
            var i = super._dimensions
            while (i-- > 0) {
                distance += Math.abs(arr[--offset] - location[i])
            }
            return distance
        }
    }

    class WeightedManhattan<T>(dims: Int) : KDTree<T>(dims) {
        private var weights: DoubleArray
        fun setWeights(newWeights: DoubleArray) {
            weights = newWeights
        }

        override fun pointRectDist(offset: Int, location: DoubleArray): Double {
            var offset = offset
            offset *= 2 * super._dimensions
            var distance = 0.0
            val array = super.nodeMinMaxBounds.array
            var i = 0
            while (i < location.size) {
                var diff = 0.0
                var bv = array[offset]
                val lv = location[i]
                if (bv > lv) diff = bv - lv else {
                    bv = array[offset + 1]
                    if (lv > bv) diff = lv - bv
                }
                distance += diff * weights[i]
                i++
                offset += 2
            }
            return distance
        }

        override fun pointDist(arr: DoubleArray, location: DoubleArray, index: Int): Double {
            var distance = 0.0
            var offset = (index + 1) * super._dimensions
            var i = super._dimensions
            while (i-- > 0) {
                distance += Math.abs(arr[--offset] - location[i]) * weights[i]
            }
            return distance
        }

        init {
            weights = DoubleArray(dims)
            for (i in 0 until dims) weights[i] = 1.0
        }
    }

    data class SearchResult<S> internal constructor(var distance: Double, var payload: S)

    private inner class Node @JvmOverloads internal constructor(pointMemory: DoubleArray = DoubleArray(_bucketSize * _dimensions)) {
        // for accessing bounding box data
// - if trees weren't so unbalanced might be better to use an implicit
// heap?
        var index: Int
        // keep track of size of subtree
        var entries = 0
        // leaf
        var pointLocations: ContiguousDoubleArrayList?
        var pointPayloads: ArrayList<T>? = ArrayList(_bucketSize)
        // stem
// Node less, more;
        var lessIndex = 0
        var moreIndex = 0
        var splitDim = 0
        var splitVal = 0.0
        fun search(searchLocation: DoubleArray, stack: IntStack) {
            if (searchLocation[splitDim] < splitVal) stack.push(moreIndex).push(lessIndex) // less will be popped
            else stack.push(lessIndex).push(moreIndex) // more will be popped
            // first
        }

        // returns number of points added to results
        fun search(searchLocation: DoubleArray, results: PrioQueue<T>): Int {
            var updated = 0
            var j = entries
            while (j-- > 0) {
                val distance = pointDist(pointLocations!!.array, searchLocation, j)
                if (results.peekPrio() > distance) {
                    updated++
                    results.addNoGrow(pointPayloads!![j], distance)
                }
            }
            return updated
        }

        fun searchBall(
            searchLocation: DoubleArray,
            radius: Double,
            results: ArrayList<T>
        ) {
            var j = entries
            while (j-- > 0) {
                val distance = pointDist(pointLocations!!.array, searchLocation, j)
                if (radius >= distance) {
                    results.add(pointPayloads!![j])
                }
            }
        }

        fun searchRect(mins: DoubleArray, maxs: DoubleArray, results: ArrayList<T>) {
            var j = entries
            while (j-- > 0) {
                if (contains(pointLocations!!.array, mins, maxs, j)) results.add(pointPayloads!![j])
            }
        }

        fun expandBounds(location: DoubleArray) {
            entries++
            var mio = index * 2 * _dimensions
            for (i in 0 until _dimensions) {
                nodeMinMaxBounds.array[mio] = Math.min(nodeMinMaxBounds.array[mio], location[i])
                mio++
                nodeMinMaxBounds.array[mio] = Math.max(nodeMinMaxBounds.array[mio], location[i])
                mio++
            }
        }

        fun add(location: DoubleArray, load: T): Int {
            pointLocations!!.add(location)
            pointPayloads!!.add(load)
            return entries
        }

        fun split() {
            var offset = index * 2 * _dimensions
            var diff = 0.0
            for (i in 0 until _dimensions) {
                val min = nodeMinMaxBounds.array[offset]
                val max = nodeMinMaxBounds.array[offset + 1]
                if (max - min > diff) {
                    var mean = 0.0
                    for (j in 0 until entries) mean += pointLocations!!.array[i + _dimensions * j]
                    mean = mean / entries
                    var varianceSum = 0.0
                    for (j in 0 until entries) varianceSum += sqr(mean - pointLocations!!.array[i + _dimensions * j])
                    if (varianceSum > diff * entries) {
                        diff = varianceSum / entries
                        splitVal = mean
                        splitDim = i
                    }
                }
                offset += 2
            }
            // kill all the nasties
            if (splitVal == Double.POSITIVE_INFINITY) splitVal =
                Double.MAX_VALUE else if (splitVal == Double.NEGATIVE_INFINITY) splitVal =
                Double.MIN_VALUE else if (splitVal == nodeMinMaxBounds.array[index * 2 * _dimensions + 2 * splitDim + 1]) splitVal =
                nodeMinMaxBounds.array[index * 2 * _dimensions + 2 * splitDim]
            val less = Node(mem_recycle) // recycle that memory!
            val more = Node()
            lessIndex = less.index
            moreIndex = more.index
            // reduce garbage by factor of _bucketSize by recycling this array
            val pointLocation = DoubleArray(_dimensions)
            for (i in 0 until entries) {
                System.arraycopy(pointLocations!!.array, i * _dimensions, pointLocation, 0, _dimensions)
                val load = pointPayloads!![i]
                if (pointLocation[splitDim] < splitVal) {
                    less.expandBounds(pointLocation)
                    less.add(pointLocation, load)
                } else {
                    more.expandBounds(pointLocation)
                    more.add(pointLocation, load)
                }
            }
            if (less.entries * more.entries == 0) { // one of them was 0, so the split was worthless. throw it away.
                _nodes -= 2 // recall that bounds memory
                nodeList.removeAt(moreIndex)
                nodeList.removeAt(lessIndex)
            } else { // we won't be needing that now, so keep it for the next split
// to reduce garbage
                mem_recycle = pointLocations!!.array
                pointLocations = null
                pointPayloads!!.clear()
                pointPayloads = null
            }
        }

        init {
            pointLocations = ContiguousDoubleArrayList(pointMemory)
            index = _nodes++
            nodeList.add(this)
            nodeMinMaxBounds.add(bounds_template)
        }
    }

    // NB! This Priority Queue keeps things with the LOWEST priority.
// If you want highest priority items kept, negate your values
    private class PrioQueue<S> internal constructor(size: Int) {
        var elements: Array<Any?>
        var priorities: DoubleArray
        private var minPrio = 0.0
        private val size: Int

        init {
            elements = arrayOfNulls(size)
            priorities = DoubleArray(size)
            Arrays.fill(priorities, Double.POSITIVE_INFINITY)
            minPrio = Double.POSITIVE_INFINITY
            this.size = size
        }


        // uses O(log(n)) comparisons and one big shift of size O(N)
// and is MUCH simpler than a heap --> faster on small sets, faster JIT
        fun addNoGrow(value: S, priority: Double) {
            val index = searchFor(priority)
            val nextIndex = index + 1
            val length = size - nextIndex
            System.arraycopy(elements, index, elements, nextIndex, length)
            System.arraycopy(priorities, index, priorities, nextIndex, length)
            elements[index] = value
            priorities[index] = priority
            minPrio = priorities[size - 1]
        }

        fun searchFor(priority: Double): Int {
            var i = size - 1
            var j = 0
            while (i >= j) {
                val index = i + j ushr 1
                if (priorities[index] < priority) j = index + 1 else i = index - 1
            }
            return j
        }

        fun peekPrio(): Double {
            return minPrio
        }

    }

    private class ContiguousDoubleArrayList internal constructor(var array: DoubleArray) {
        var size = 0

        internal constructor(size: Int) : this(DoubleArray(size)) {}

        fun add(da: DoubleArray): ContiguousDoubleArrayList {
            if (size + da.size > array.size) array = Arrays.copyOf(array, (array.size + da.size) * 2)
            System.arraycopy(da, 0, array, size, da.size)
            size += da.size
            return this
        }

    }

    private class IntStack internal constructor(var array: IntArray) {
        var size = 0

        @JvmOverloads
        internal constructor(size: Int = 64) : this(IntArray(size)) {
        }

        fun push(i: Int): IntStack {
            if (size >= array.size) array = Arrays.copyOf(array, (array.size + 1) * 2)
            array[size++] = i
            return this
        }

        fun pop(): Int {
            return array[--size]
        }

        fun size(): Int {
            return size
        }

    }

    companion object {
        // use a big bucketSize so that we have less node bounds (for more cache
// hits) and better splits
// if you have lots of dimensions this should be big, and if you have few small
        private const val _bucketSize = 50

        fun sqr(d: Double): Double {
            return d * d
        }
    }

    init {
        // initialise this big so that it ends up in 'old' memory
        nodeMinMaxBounds = ContiguousDoubleArrayList(512 * 1024 / 8 + 2 * _dimensions)
        mem_recycle = DoubleArray(_bucketSize * _dimensions)
        bounds_template = DoubleArray(2 * _dimensions)
        Arrays.fill(bounds_template, Double.NEGATIVE_INFINITY)
        var i = 0
        val max = 2 * _dimensions
        while (i < max) {
            bounds_template[i] = Double.POSITIVE_INFINITY
            i += 2
        }
        // and.... start!
        root = Node()
    }
}