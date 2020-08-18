package hermiteForm

import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import org.apache.commons.math3.exception.OutOfRangeException
import org.apache.commons.math3.util.ArithmeticUtils
import org.apache.commons.math3.util.ArithmeticUtils.gcd
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashMap
import kotlin.math.*

open class SparseColIntMatrix: ArrayList<SparseColIntMatrix.SparseIntColumn> {
    var nRows: Int
    val nCols: Int
        get() = this.size

    val columns: MutableList<SparseIntColumn>
        get() = this

    constructor(nRows: Int, nCols: Int): super(nCols) {
        for(col in 1..nCols) add(SparseIntColumn())
        this.nRows = nRows
    }

    constructor(copy: SparseColIntMatrix): this(copy, copy.nRows)

    constructor(copy: Collection<SparseIntColumn>, nRows: Int): super(copy.size) {
        for(col in copy) add(SparseIntColumn(col))
        this.nRows = nRows
    }

    operator fun set(row: Int, col: Int, element: Int) {
        this[col][row] = element
    }

    operator fun get(row: Int, col: Int): Int = this[col][row]

    fun swapCols(col1: Int, col2: Int) {
        val tmp = this[col1]
        this[col1] = this[col2]
        this[col2] = tmp
    }

    fun setRow(row: Int, rowEntries: Collection<Map.Entry<Int,Int>>) {
        if(row >= nRows) throw(OutOfRangeException(row, 0, nRows-1))
        for(entry in rowEntries) {
            this[entry.key][row] = entry.value
        }
    }

    fun addRow(rowEntries: Collection<Map.Entry<Int,Int>>) {
        for(entry in rowEntries) {
            this[entry.key][nRows] = entry.value
        }
        nRows++
    }

    fun addColLeft(column: SparseIntColumn) {
        this.add(0,column)
    }

    fun addColRight(column: SparseIntColumn) {
        this.add(column)
    }

    fun transpose(): SparseColIntMatrix {
        val transpose = SparseColIntMatrix(nCols, nRows)
        for(thisCol in 0 until nCols) {
            for(entry in this[thisCol].entries) {
                transpose[thisCol, entry.key] = entry.value
            }
        }
        return transpose
    }


    fun nonZeroElementCount(): Int = sumBy { it.entries.size }


    // solves LX=Y where L is the lower triangular
    // part of this matrix
    fun solveLowerTriangular(Y: SparseIntColumn): SparseIntColumn {
        val X = SparseIntColumn(Y)
        for(i in 0 until min(nCols, nRows)) {
            val Xi = X[i]/this[i,i]
//            assert(Xi*this[i,i] == X[i])
            X[i] = Xi
            for(entry in this[i].entries) {
                if(entry.key > i) X[entry.key] -= Xi*entry.value
            }
        }
        return X
    }


    operator fun times(X: SparseIntColumn): SparseIntColumn {
        val result = SparseIntColumn()
        for(entry in X.entries) {
//            assert(entry.value != 0)
            result.weightedPlusAssign(this[entry.key], entry.value)
        }
        return result
    }

    operator fun times(X: IntArray): SparseIntColumn {
        val result = SparseIntColumn()
        for(i in X.indices) {
            if(X[i] != 0) result.weightedPlusAssign(this[i], X[i])
        }
        return result
    }


    operator fun times(X: List<Int>): SparseIntColumn {
        val result = SparseIntColumn()
        for(i in X.indices) {
            if(X[i] != 0) result.weightedPlusAssign(this[i], X[i])
        }
        return result
    }


    fun copy(): SparseColIntMatrix = SparseColIntMatrix(this)


    fun hermiteDecomposition(): SparseColIntMatrix {
        val U = Identity(this.nCols)
        for(eqnIndex in 0 until this.nRows) {
            var indexCol = this[eqnIndex]
            var UindexCol = U[eqnIndex]
            var indexColMagnitude = indexCol[eqnIndex]
            for(colIndex in eqnIndex+1 until this.nCols) {
                val swapCol = this[colIndex]
                val UswapCol = U[colIndex]
                val swapColMagnitude = swapCol[eqnIndex]
                if(swapColMagnitude != 0) {
                    if(indexColMagnitude != 0) {
                        val g = ExtendedEuclid(indexColMagnitude, swapColMagnitude)
                        val indexMultiplier = swapColMagnitude/g.gcd
                        val swapColMultiplier = indexColMagnitude/g.gcd
                        hermiteColumnSwap(indexCol, swapCol, g.y, swapColMultiplier, indexMultiplier)
                        hermiteColumnSwap(UindexCol, UswapCol, g.y, swapColMultiplier, indexMultiplier)
                        indexColMagnitude = indexCol[eqnIndex]
                    } else {
                        this.swapCols(eqnIndex, colIndex)
                        indexCol = this[eqnIndex]
                        U.swapCols(eqnIndex, colIndex)
                        UindexCol = U[eqnIndex]
                        indexColMagnitude = swapColMagnitude
                    }
                }
            }
        }
        return U
    }

    fun upperTriangularise() {
        var col = nCols-1
        var row = nRows-1
        while(row >= 0 && col > 0) {
            var indexCol = this[col]
            var indexColMagnitude: Int
            do {
                indexColMagnitude = indexCol[row]
                for (colIndex in col - 1 downTo 0) {
                    val swapCol = this[colIndex]
                    val swapColMagnitude = swapCol[row]
                    if (swapColMagnitude != 0) {
                        if (indexColMagnitude != 0) {
                            val g = ExtendedEuclid(indexColMagnitude, swapColMagnitude)
                            hermiteColumnSwap(indexCol, swapCol, g.y, indexColMagnitude/g.gcd, swapColMagnitude/g.gcd)
//                            val g = gcd(indexColMagnitude, swapColMagnitude)
//                            swapCol *= -indexColMagnitude / g
//                            swapCol.weightedPlusAssign(indexCol, swapColMagnitude / g)
                        } else {
                            this.swapCols(col, colIndex)
                            indexCol = this[col]
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
            val lastCol = this[col]
            var lastColMagnitude: Int
            do {
                lastColMagnitude = lastCol[--row]
            } while (row >= 0 && lastColMagnitude == 0)
            if (lastColMagnitude < 0) lastCol *= -1
        }
    }


    fun upperTriangularThin() {
        for(col in nCols-2 downTo 0) {
            val indexCol = this[col]
            indexCol.keys.max()?.also { row ->
                val indexColMagnitude = this[row,col]
                for(j in col+1 until nCols) {
                    val swapCol = this[j]
                    val swapColMagnitude = swapCol[row]
                    if(swapColMagnitude != 0) {
                        val g = gcd(indexColMagnitude, swapColMagnitude)
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
        val activeCols = LinkedList(this)
        for(row in 0 until nRows) {
            var indexCol: SparseIntColumn? = null
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
            val reducingCol = this[reducingi]
            val reducingColNormSq = reducingCol.normSqL2()
            for(i in reducingi+1 until nCols) {
                val colToReduce = this[i]
                val n = (reducingCol.dotProd(colToReduce).toDouble()/reducingColNormSq).roundToInt()
//                println("n is $n")
                colToReduce.weightedPlusAssign(reducingCol, -n)
            }
        }
    }

    fun minimiseAfterDoubleTriangularisation() {
        // first calculate ranges
        val ranges = Array(nCols) { i ->
            IntRange(this[i].keys.min()?:0, this[i].keys.max()?:0)
        }
        for(reducingi in 0 until nCols-1) {
            val reducingCol = this[reducingi]
            val reducingColNormSq = reducingCol.normSqL2()
            for(i in reducingi+1 until nCols) {
                if(ranges[i].first < ranges[reducingi].first && ranges[i].last > ranges[reducingi].last) { // TODO: reducce this to isSubsetOf relation
                    val colToReduce = this[i]
                    val n = (reducingCol.dotProd(colToReduce).toDouble()/reducingColNormSq).roundToInt()
                    colToReduce.weightedPlusAssign(reducingCol, -n)
                }
            }
        }
    }

    // Minimise CX
    // Subject to:
    // MX >= B
    // and
    // X >= 0 and X is integer
    // where M is this matrix
    // returns X
    //
    fun IPsolve(B: SparseIntColumn, C: List<Double> = DoubleArray(nCols) {0.0}.asList()): IntArray {
        val solver = MPSolver("SparseSolver", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val X = solver.makeIntVarArray(nCols, 0.0, Double.POSITIVE_INFINITY)
        val constraints = Array<MPConstraint>(nRows) { solver.makeConstraint() }
        for(col in 0 until nCols) {
            for(coefficient in this[col].entries) {
                constraints[coefficient.key].setCoefficient(X[col], coefficient.value.toDouble())
            }
        }
        for(i in 0 until nRows) {
            constraints[i].setBounds(B[i].toDouble(), Double.POSITIVE_INFINITY)
        }
        val objective = solver.objective()
        for(i in C.indices) {
            objective.setCoefficient(X[i], C[i])
        }
        solver.objective().setMinimization()
        val solveState = solver.solve()
        return if (solveState == MPSolver.ResultStatus.OPTIMAL)
            IntArray(X.size) { i -> X[i].solutionValue().toInt() }
        else
            throw(RuntimeException(
                when (solveState) {
                    MPSolver.ResultStatus.INFEASIBLE -> "Infeasible"
                    MPSolver.ResultStatus.UNBOUNDED -> "Unbounded"
                    else -> "Solve Error"
                }
            ))
    }


    // sets
    // c1' = X*c1 + Y*c2
    // c2' = B*c1 - A*c2
    // where X is an integer, A != 0 and XA + YB = 1
    // which is equivalent to
    // c2 = B*c1 - A*c2
    // c1 = (c1 - Y*c2)/A
    private fun hermiteColumnSwap(c1: SparseIntColumn, c2: SparseIntColumn, Y: Int, A: Int, B: Int) {
        c2 *= -A
        c2.weightedPlusAssign(c1,B)
// works better without reweighting the index column
        c1.weightedPlusAssign(c2,-Y)
        c1 /= A
    }



    override fun toString(): String {
        val out = StringBuilder()
        for(row in 0 until nRows) {
            for (col in this) {
                val v = col[row]
                out.append(v)
                out.append(' ')
            }
            out.appendln()
        }
        return out.toString()
    }

    fun toSparsityString(): String {
        val out = StringBuilder()
        for(row in 0 until nRows) {
            for (col in this) {
                val v = col[row]
                out.append(when(v.sign) {
                    1 -> '+'
                    -1 -> '-'
                    else ->  '.'
                })
            }
            out.appendln()
        }
        return out.toString()
    }


    // map from column index to non-zero entry value
    class SparseIntColumn(val data: HashMap<Int,Int> = HashMap(2)) {
        val entries: MutableSet<MutableMap.MutableEntry<Int, Int>>
            get() = data.entries

        val keys: MutableSet<Int>
            get() = data.keys

        val values: MutableCollection<Int>
            get() = data.values

        constructor(copy: SparseIntColumn): this(HashMap(copy.data))

        constructor(denseCol: IntArray): this() {
            for(i in 0 until denseCol.size) {
                val v = denseCol[i]
                if(v != 0) data[i] = v
            }
        }

        constructor(vararg entries: Pair<Int,Int>): this(HashMap(mapOf(*entries)))

        operator fun get(key: Int): Int {
            return data.getOrDefault(key, 0)
        }

        operator fun set(key: Int, value: Int) {
            if(value != 0) data[key] = value else data.remove(key)
        }

        operator fun plus(other: SparseIntColumn): SparseIntColumn {
            val sum = SparseIntColumn(this)
            for(entry in other.entries) {
                sum.data.merge(entry.key, entry.value) {a,b ->
                    val result = a+b
                    if(result !=0) result else null
                }
            }
            return sum
        }

        operator fun minus(other: SparseIntColumn): SparseIntColumn {
            val sum = SparseIntColumn(this)
            sum -= other
            return sum
        }

        operator fun timesAssign(multiplier: Int) {
            if(multiplier == 0) data.clear()
            for(entry in entries) {
                entry.setValue(entry.value * multiplier)
            }
        }

        operator fun divAssign(denominator: Int) {
            for(entry in entries) {
                entry.setValue(entry.value/denominator)
            }
        }

        operator fun minusAssign(otherCol: SparseIntColumn) {
            for(entry in otherCol.entries) {
                data.merge(entry.key, -entry.value) { a, b ->
                    val result = a + b
                    if(result != 0) result else null
                }
            }
        }

        // this += weight*otherCol
        fun weightedPlusAssign(otherCol: SparseIntColumn, weight: Int) {
            if(weight == 0) return
            for(entry in otherCol.entries) {
                data.merge(entry.key, entry.value*weight) { a, b ->
                    val result = a + b
                    if(result != 0) result else null
                }
            }
        }


        fun toIntArray(size: Int): IntArray {
            val array = IntArray(size)
            for(entry in this.data) {
                array[entry.key] = entry.value
            }
            return array
        }


        fun isPositive(): Boolean {
            return data.values.fold(true) { isPositive, v -> isPositive && (v >= 0) }
        }

        fun normSqL2(): Int = this.values.sumBy { it*it }

        fun dotProd(other: SparseIntColumn): Int {
            val columns =
                if(this.entries.size < other.entries.size) Pair(this, other) else Pair(other,this)
            return columns.first.entries.sumBy { it.value * columns.second[it.key] }
        }

        fun copy() = SparseIntColumn(this)

        fun toTrajectory(dynamics: SparseColIntMatrix, agentRows: Int, agentCols: Int): String {
            val trajectory = Array(agentCols) { Array(agentRows) { 0 } }
            for(act in entries) {
                val particleVal = if(act.value > 0) 1 else 2
                for(actPosition in dynamics[act.key].keys) {
                    val x = actPosition.rem(agentCols)
                    val y = actPosition.div(agentCols).rem(agentRows)
                    trajectory[x][y] = trajectory[x][y] or particleVal
                }
            }
            val particleChars = arrayOf(' ','+','-','*')
            val trajectoryString = StringBuilder()
            for(row in agentRows-1 downTo 0) {
                for(col in 0 until agentCols) {
                    trajectoryString.append(particleChars[trajectory[row][col]])
                }
                trajectoryString.appendln()
            }
            return trajectoryString.toString()
        }

        override fun toString(): String {
            return data.toString()
        }

        operator fun unaryMinus(): SparseColIntMatrix.SparseIntColumn {
            val negation = SparseIntColumn()
            negation.weightedPlusAssign(this, -1)
            return negation
        }
    }

    class ExtendedEuclid {
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
//            x = last_x
            y = last_y
            gcd = a
//            assert(A.rem(gcd) == 0)
//            assert(B.rem(gcd) == 0)
//            assert(A*last_x + B*last_y == gcd)
        }
    }

    companion object {
        fun Identity(size: Int): SparseColIntMatrix {
            val I = SparseColIntMatrix(size,size)
            for(i in 0 until size) {
                I[i,i] = 1
            }
            return I
        }
    }
}