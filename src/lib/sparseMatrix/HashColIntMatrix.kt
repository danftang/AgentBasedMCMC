package lib.sparseMatrix

import org.apache.commons.math3.exception.OutOfRangeException
import org.apache.commons.math3.util.ArithmeticUtils.gcd
import java.util.*
import kotlin.collections.ArrayList
import kotlin.math.*

open class HashColIntMatrix: SparseColIntMatrix {
    override var nRows: Int
    override val nCols: Int
        get() = _columns.size
    override val columns: List<SparseIntVector>
        get() = _columns

    private val _columns: ArrayList<HashIntVector>


    constructor(nRows: Int, nCols: Int) {
        _columns = ArrayList(nCols)
        for(col in 1..nCols) _columns.add(HashIntVector())
        this.nRows = nRows
    }


    constructor(copy: SparseIntMatrix): this(copy.nRows, copy.nCols) {
        for(entry in copy.entries) {
            _columns[entry.col][entry.row] = entry.value
        }
    }

    override operator fun set(row: Int, col: Int, value: Int) {
        if(value == 0) {
            _columns[col].data.remove(row)
            return
        }
        _columns[col][row] = value
    }

    override operator fun get(row: Int, col: Int): Int = columns[col][row]

    override fun plusAssign(row: Int, col: Int, addition: Int) {
        _columns[col].data.merge(row, addition, Int::plus)
    }

    override fun clearColumn(j: Int) {
        _columns[j].data.clear()
    }

    override fun removeColumns(colsToRemove: Iterable<Int>) {
        _columns.removeAll(colsToRemove.map { _columns[it] })
    }

    override fun replaceNonZeroElementsInCol(col: Int, map: (row: Int, value: Int) -> Int) {
        _columns[col].data.entries.forEach { entry ->
            entry.setValue(map(entry.key, entry.value))
        }
    }

    override fun resize(nRows: Int, nCols: Int) {
        if(nCols > columns.size) {
            repeat(nCols - columns.size) { _columns.add(HashIntVector()) }
        } else if(nCols < columns.size) {
            repeat(columns.size - nCols) {
                _columns.removeAt(_columns.lastIndex)
            }
        }
        if(nRows < this.nRows) {
            for(col in _columns) {
                col.data.keys.removeIf { it >= nRows }
            }
        }
        this.nRows = nRows
    }

    override fun swapCols(col1: Int, col2: Int) {
        val tmp = _columns[col1]
        _columns[col1] = _columns[col2]
        _columns[col2] = tmp
    }

    fun setRow(row: Int, rowEntries: Collection<Map.Entry<Int,Int>>) {
        if(row >= nRows) throw(OutOfRangeException(row, 0, nRows-1))
        for(entry in rowEntries) {
            this[row,entry.key] = entry.value
        }
    }

//    fun addRow(rowEntries: Collection<Map.Entry<Int,Int>>) {
//        for(entry in rowEntries) {
//            this[entry.key][nRows] = entry.value
//        }
//        nRows++
//    }
//
//    fun addColLeft(column: SparseIntColumn) {
//        this.add(0,column)
//    }
//
//    fun addColRight(column: SparseIntColumn) {
//        this.add(column)
//    }

    fun transpose(): HashColIntMatrix {
        val transpose = HashColIntMatrix(nCols, nRows)
        for(thisCol in 0 until nCols) {
            for(entry in columns[thisCol]) {
                transpose[thisCol, entry.key] = entry.value
            }
        }
        return transpose
    }



    fun nonZeroElementCount(): Int = columns.sumBy { it.sparseSize }


    operator fun times(M: SparseIntMatrix): HashColIntMatrix {
        assert(nCols == M.nRows)
        val product = HashColIntMatrix(nRows, M.nCols)
        for(otherEntry in M.entries) {
            for (entry in columns[otherEntry.row]) {
                product[entry.key, otherEntry.col] += entry.value * otherEntry.value
            }
        }
        return product
    }

    fun copy(): HashColIntMatrix = HashColIntMatrix(this)

    fun subMatrixView(fromColumn: Int, untilColumn: Int): HashColIntMatrix {
        val view = HashColIntMatrix(nRows, untilColumn - fromColumn)
        view._columns.ensureCapacity(untilColumn - fromColumn)
        for(j in fromColumn until untilColumn) {
            view._columns.add(_columns[j])
        }
        return view
    }

    // Number of non-zero elements / total number of elements
    fun sparsityRatio(): Double {
        return columns.sumBy { it.sparseSize }.toDouble() / (nRows.toDouble()*nCols)
    }

//    fun hermiteDecomposition(): HashColIntMatrix {
//        val U = identity(this.nCols)
//        for(eqnIndex in 0 until this.nRows) {
////            println("Calculating row $eqnIndex")
//            var indexCol = _columns[eqnIndex]
//            var UindexCol = U._columns[eqnIndex]
//            var indexColMagnitude = indexCol[eqnIndex]
//            for(colIndex in this.nCols-1 downTo eqnIndex+1) {
//                val swapCol = _columns[colIndex]
//                val UswapCol = U._columns[colIndex]
//                val swapColMagnitude = swapCol[eqnIndex]
//                if(swapColMagnitude != 0) {
//                    if(indexColMagnitude != 0) {
//                        val g = ExtendedEuclid(indexColMagnitude, swapColMagnitude)
//                        val indexMultiplier = swapColMagnitude/g.gcd
//                        val swapColMultiplier = indexColMagnitude/g.gcd
//                        hermiteColumnSwap(indexCol, swapCol, g.y, swapColMultiplier, indexMultiplier)
//                        hermiteColumnSwap(UindexCol, UswapCol, g.y, swapColMultiplier, indexMultiplier)
//                        indexColMagnitude = indexCol[eqnIndex]
//                    } else {
//                        this.swapCols(eqnIndex, colIndex)
//                        indexCol = _columns[eqnIndex]
//                        U.swapCols(eqnIndex, colIndex)
//                        UindexCol = U._columns[eqnIndex]
//                        indexColMagnitude = swapColMagnitude
//                    }
//                }
//            }
//            // now reduce to the left
////            for(colIndex in 0 until eqnIndex) {
////                val swapColMagnitude = this[eqnIndex,colIndex]
////                if(swapColMagnitude != 0) {
////                    val reductionWeight = -swapColMagnitude / indexColMagnitude
////                    this[colIndex].weightedPlusAssign(indexCol, reductionWeight)
////                    U[colIndex].weightedPlusAssign(UindexCol, reductionWeight)
////                }
////            }
//        }
//        return U
//    }


    fun upperTriangularise() {
        var col = nCols-1
        var row = nRows-1
        while(row >= 0 && col > 0) {
//            println("upper triangularising row $row")
            var indexCol = _columns[col]
            var indexColMagnitude: Int
            do {
                indexColMagnitude = indexCol[row]
                for (colIndex in col - 1 downTo 0) {
                    val swapCol = _columns[colIndex]
                    val swapColMagnitude = swapCol[row]
                    if (swapColMagnitude != 0) {
                        if (indexColMagnitude != 0) {
                            val g = ExtendedEuclid(indexColMagnitude, swapColMagnitude)
                            hermiteColumnSwap(indexCol, swapCol, g.y, indexColMagnitude/g.gcd, swapColMagnitude/g.gcd)
                            indexColMagnitude = indexCol[row]
//                            val g = gcd(indexColMagnitude, swapColMagnitude)
//                            swapCol *= -indexColMagnitude / g
//                            swapCol.weightedPlusAssign(indexCol, swapColMagnitude / g)
                        } else {
                            this.swapCols(col, colIndex)
                            indexCol = _columns[col]
                            indexColMagnitude = swapColMagnitude
                        }
                    }
                }
                --row
            } while(indexColMagnitude == 0 && row >= 0)
            if(indexColMagnitude < 0) indexCol *= -1 // ensure leading element is positive
            --col
        }
        if(row > 0 && col >=0) {
            // ensure leading element of last column is positive
            val lastCol = _columns[col]
            var lastColMagnitude: Int
            do {
                lastColMagnitude = lastCol[--row]
            } while (row >= 0 && lastColMagnitude == 0)
            if (lastColMagnitude < 0) lastCol *= -1
        }
    }


    fun upperTriangularThin() {
        for(col in nCols-2 downTo 0) {
//            println("Thinning with col $col")
            val indexCol = columns[col]
            indexCol.keys.max()?.also { row ->
                val indexColMagnitude = this[row,col]
                for(j in col+1 until nCols) {
                    val swapCol = _columns[j]
                    val swapColMagnitude = swapCol[row]
                    if(swapColMagnitude != 0) {
                        val g = if(indexColMagnitude == 1) 1 else gcd(indexColMagnitude, swapColMagnitude)
                        if(indexColMagnitude > 0) {
                            swapCol *= indexColMagnitude / g
                            swapCol.weightedPlusAssign(indexCol, -swapColMagnitude / g)
                        } else {
                            swapCol *= -indexColMagnitude / g
                            swapCol.weightedPlusAssign(indexCol, swapColMagnitude / g)
                        }
                    }
                }
            }
        }
    }

    fun lowerTriangulariseWithoutSwap() {
        val activeCols = LinkedList(_columns)
        for(row in 0 until nRows) {
            var indexCol: SparseIntVector? = null
            var indexMagnitude = 0
            val colIt = activeCols.iterator()
            while(colIt.hasNext()) {
                val col = colIt.next()
                if(col[row] != 0) {
                    if(indexCol == null) { // must be first occurence of non-zero in this row
                        indexCol = col
                        indexMagnitude = indexCol[row]
                        colIt.remove()
                    } else {
                        val swapColMagnitude = col[row]
                        val g = gcd(indexMagnitude, swapColMagnitude)
                        col *= -indexMagnitude / g
                        col.weightedPlusAssign(indexCol, swapColMagnitude / g)
                    }
                }
            }
        }
    }

    fun minimiseL2NormLeftToRight() {
        for(reducingi in 0 until nCols-1) {
            val reducingCol = columns[reducingi]
            val reducingColNormSq = reducingCol.normSqL2()
            for(i in reducingi+1 until nCols) {
                val colToReduce = _columns[i]
                val n = (reducingCol.dotProd(colToReduce).toDouble()/reducingColNormSq).roundToInt()
                colToReduce.weightedPlusAssign(reducingCol, -n)
            }
        }
    }

    fun minimiseAfterDoubleTriangularisation() {
        // first calculate ranges
        val ranges = Array(nCols) { i ->
            IntRange(columns[i].keys.min()?:0, columns[i].keys.max()?:0)
        }
        for(reducingi in 0 until nCols-1) {
            val reducingCol = columns[reducingi]
            val reducingColNormSq = reducingCol.normSqL2()
            for(i in reducingi+1 until nCols) {
                if(ranges[i].first < ranges[reducingi].first && ranges[i].last > ranges[reducingi].last) { // TODO: reducce this to isSubsetOf relation
                    val colToReduce = _columns[i]
                    val n = (reducingCol.dotProd(colToReduce).toDouble()/reducingColNormSq).roundToInt()
                    colToReduce.weightedPlusAssign(reducingCol, -n)
                }
            }
        }
    }


    // sets
    // c1' = X*c1 + Y*c2
    // c2' = B*c1 - A*c2
    // where X is an integer, A != 0 and XA + YB = 1
    // which is equivalent to
    // c2' = B*c1 - A*c2
    // c1' = (c1 - Y*c2')/A
    private fun hermiteColumnSwap(c1: HashIntVector, c2: HashIntVector, Y: Int, A: Int, B: Int) {
        c2 *= -A
        c2.weightedPlusAssign(c1,B)
// works better without reweighting the index column
        c1.weightedPlusAssign(c2,-Y)
        c1 /= A
    }


    override fun toString(): String {
        val out = StringBuilder()
        for(row in 0 until nRows) {
            for (col in columns) {
                val v = col[row]
                if(v in 0..9) out.append(' ')
                out.append(v)
                out.append(' ')
            }
            out.append('\n')
        }
        return out.toString()
    }


    fun checkNoZeroEntries(): Boolean {
        return columns.all { col -> col.all { it.value != 0 } }
    }

    fun valueRange(): IntRange {
        var max = Int.MIN_VALUE
        var min = Int.MAX_VALUE
        _columns.forEach { col ->
            val colRange = col.valueRange()
            if(colRange.start < min) min = colRange.start
            if(colRange.endInclusive > max) max = colRange.endInclusive
        }
        return IntRange(min, max)
    }

    // map from column index to non-zero entry value
//    class SparseIntColumn(val data: HashMap<Int,Int> = HashMap(4)): SparseIntVector {
//        val entries: MutableSet<MutableMap.MutableEntry<Int, Int>>
//            get() = data.entries
//
//        override val keys: MutableSet<Int>
//            get() = data.keys
//
//        override val values: MutableCollection<Int>
//            get() = data.values
//
//        override val sparseSize: Int
//            get() = entries.size
//
//        constructor(copy: SparseIntColumn): this(HashMap(copy.data))
//
//        constructor(denseCol: IntArray): this() {
//            for(i in 0 until denseCol.size) {
//                val v = denseCol[i]
//                if(v != 0) data[i] = v
//            }
//        }
//
//        constructor(vararg entries: Pair<Int,Int>): this(HashMap(mapOf(*entries)))
//
//        override operator fun get(key: Int): Int {
//            return data.getOrDefault(key, 0)
//        }
//
//        operator fun set(key: Int, value: Int) {
//            if(value != 0) data[key] = value else data.remove(key)
//        }
//
//        operator fun plus(other: SparseIntColumn): SparseIntColumn {
//            val sum = SparseIntColumn(this)
//            for(entry in other.entries) {
//                sum.data.merge(entry.key, entry.value) {a,b ->
//                    val result = a+b
//                    if(result !=0) result else null
//                }
//            }
//            return sum
//        }
//
//        operator fun plus(other: IntArray): SparseIntColumn {
//            val sum = SparseIntColumn(this)
//            for(i in other.indices) {
//                if(other[i] != 0) {
//                    sum.data.merge(i, other[i]) { a, b ->
//                        val result = a + b
//                        if (result != 0) result else null
//                    }
//                }
//            }
//            return sum
//        }
//
//
//        operator fun minus(other: SparseIntColumn): SparseIntColumn {
//            val sum = SparseIntColumn(this)
//            sum -= other
//            return sum
//        }
//
//        operator fun timesAssign(multiplier: Int) {
//            if(multiplier == 0) data.clear()
//            if(multiplier == 1) return
//            for(entry in entries) {
//                entry.setValue(entry.value * multiplier)
//            }
//        }
//
//        operator fun divAssign(denominator: Int) {
//            if(denominator == 1) return
//            for(entry in entries) {
//                entry.setValue(entry.value/denominator)
//            }
//        }
//
//        operator fun minusAssign(otherCol: SparseIntColumn) {
//            for(entry in otherCol.entries) {
//                data.merge(entry.key, -entry.value) { a, b ->
//                    val result = a + b
//                    if(result != 0) result else null
//                }
//            }
//        }
//
//        override fun equals(other: Any?): Boolean {
//            return if(other is SparseIntColumn) {
//                data == other.data
//            } else false
//        }
//
//        override fun iterator(): Iterator<Map.Entry<Int, Int>> {
//            return data.iterator()
//        }
//
//        // this += weight*otherCol
//        fun weightedPlusAssign(otherCol: SparseIntColumn, weight: Int) {
//            if(weight == 0) return
//            for(entry in otherCol.entries) {
//                data.merge(entry.key, entry.value*weight) { a, b ->
//                    val result = a + b
//                    if(result != 0) result else null
//                }
//            }
//        }
//
//
//        fun dotProd(other: SparseIntColumn): Int {
//            val columns =
//                if(this.entries.size < other.entries.size) Pair(this, other) else Pair(other,this)
//            return columns.first.entries.sumBy { it.value * columns.second[it.key] }
//        }
//
//        fun copy() = SparseIntColumn(this)
//
//        fun toTrajectory(dynamics: SparseColIntMatrix, agentRows: Int, agentCols: Int): String {
//            val trajectory = Array(agentCols) { Array(agentRows) { 0 } }
//            for(act in entries) {
//                val particleVal = if(act.value > 0) 1 else 2
//                for(actPosition in dynamics.columns[act.key].keys) {
//                    val x = actPosition.rem(agentCols)
//                    val y = actPosition.div(agentCols).rem(agentRows)
//                    trajectory[x][y] = trajectory[x][y] or particleVal
//                }
//            }
//            val particleChars = arrayOf(' ','+','-','*')
//            val trajectoryString = StringBuilder()
//            for(row in agentRows-1 downTo 0) {
//                for(col in 0 until agentCols) {
//                    trajectoryString.append(particleChars[trajectory[row][col]])
//                }
//                trajectoryString.appendln()
//            }
//            return trajectoryString.toString()
//        }
//
//        override fun toString(): String {
//            return data.toString()
//        }
//
//        operator fun unaryMinus(): SparseColIntMatrix.SparseIntColumn {
//            val negation = SparseIntColumn()
//            negation.weightedPlusAssign(this, -1)
//            return negation
//        }
//
//        fun checkNoZeroEntries(): Boolean {
//            return values.fold(true){acc, v -> acc && v != 0}
//        }
//
//    }

    class ExtendedEuclid {
        val x: Int
        val y: Int
        val gcd: Int

        constructor(A: Int, B: Int) {
            var a = A
            var b = B
            var xp = 0
            var yp = 1
            var last_x = 1
            var last_y = 0
            var temp: Int
            while (b != 0) {
                val q = a / b
                val r = a % b
                a = b
                b = r
                temp = xp
                xp = last_x - q * xp
                last_x = temp
                temp = yp
                yp = last_y - q * yp
                last_y = temp
            }
            x = last_x
            y = last_y
            gcd = a
//            assert(A.rem(gcd) == 0)
//            assert(B.rem(gcd) == 0)
//            assert(A*last_x + B*last_y == gcd)
        }
    }

    companion object {
        fun identity(size: Int): HashColIntMatrix {
            val I = HashColIntMatrix(size,size)
            for(i in 0 until size) {
                I[i,i] = 1
            }
            return I
        }
    }
}