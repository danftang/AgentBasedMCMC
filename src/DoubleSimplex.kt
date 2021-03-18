import org.apache.commons.math3.util.OpenIntToDoubleHashMap
import java.lang.RuntimeException
import kotlin.math.absoluteValue

class DoubleSimplex(constraints: List<MutableConstraint<Number>>, objective: Map<Int,Number>, initialSolution: Map<Int,Number>) {

    val epsilon = 1e-6
    val nConstraints: Int = constraints.size

    val rows = ArrayList<OpenIntToDoubleHashMap>(nConstraints + 1)

    var firstSlackColumn: Int = Integer.max(constraints.numVars(),
        objective.keys
            .max()
            ?.let { it + 1 }
            ?: 0
    )

    val nVariables: Int = firstSlackColumn + constraints.count { it.relation != "==" }

    val basicColsByRow = IntArray(nConstraints) { -1 }
    val basicRowsByCol = IntArray(nVariables) { -1 }


    val M: Matrix = object: Matrix {
        override fun get(i: Int, j: Int): Double {
            return rows[i][j]
        }

        override fun set(i: Int, j: Int, value: Double) {
            if(value.absoluteValue > epsilon) rows[i].put(j, value) else rows[i].remove(j)
        }
    }

    val B = object: Vector {
        override fun get(i: Int) = rows[i][nVariables]

        override fun set(i: Int, value: Double) {
            if(value.absoluteValue > epsilon) rows[i].put(nVariables, value) else rows[i].remove(nVariables)
        }
    }

    val objective = object: Vector {
        override fun get(j: Int) = rows[nConstraints][j]

        override fun set(j: Int, value: Double) {
            if(value.absoluteValue > epsilon) rows[nConstraints].put(j, value) else rows[nConstraints].remove(j)
        }
    }



    init {
        println("Transferring constraints to simplex tableau...")
        var nextSlackVar = firstSlackColumn
        val initialNonZeroColumns = ArrayList(initialSolution.keys)
        rows.ensureCapacity(nConstraints + 1)
        repeat(nConstraints + 1) { rows.add(OpenIntToDoubleHashMap(0.0)) }
        constraints.forEachIndexed { i, constraint ->
            val slackness = constraint.slackness(initialSolution.mapValues { it.value.toDouble() })
            if(slackness < 0.0 || (constraint.relation == "==" && slackness != 0.0)) println("WARNING: initial solution is not feasible: Slackness $slackness")
            if(slackness != 0.0) initialNonZeroColumns.add(nextSlackVar)
            if (constraint.constant.toDouble() >= 0.0) {
                constraint.coefficients.forEach { (j, x) -> M[i, j] = x.toDouble() }
                B[i] = constraint.constant.toDouble()
                when (constraint.relation) {
                    "<=" -> { basicColsByRow[i] = nextSlackVar; M[i, nextSlackVar++] = 1.0 }
                    ">=" -> { M[i, nextSlackVar++] = -1.0 }
                }
            } else {
                constraint.coefficients.forEach { (j, x) -> M[i, j] = -x.toDouble() }
                B[i] = -constraint.constant.toDouble()
                when (constraint.relation) {
                    "<=" -> { M[i, nextSlackVar++] = -1.0 }
                    ">=" -> { basicColsByRow[i] = nextSlackVar; M[i, nextSlackVar++] = 1.0 }
                }
            }
        }
        objective.forEach { (j,x) -> this.objective[j] = x.toDouble() }

        println("Initial size = ${rows.sumBy {it.size()}} ")
        println("Pivoting in initial solution")
        initialPivot(initialNonZeroColumns)
        println("Done")
    }


    fun X(includeSlacks: Boolean=false): Map<Int,Double> {
        val x = HashMap<Int,Double>()
        for (i in basicColsByRow.indices) {
            val basicCol = basicColsByRow[i]
            val Bi = B[i]
            if((includeSlacks || basicCol < firstSlackColumn) && Bi != 0.0) x[basicCol] = B[i] / M[i, basicCol]
        }
        return x
    }


    // pivots in all columns in the supplied list
    // then greedily pivots in any remaining rows
    fun initialPivot(columnsToInclude: List<Int>) {
        val canPivotRow = Array(basicColsByRow.size) { true }
        var pivotRow: Int
        for(j in columnsToInclude) {
            if(!isBasicColumn(j)) {
                pivotRow = 0
                while(pivotRow < nConstraints && (!canPivotRow[pivotRow] || M[pivotRow,j] == 0.0)) ++pivotRow
                if(pivotRow == nConstraints) throw(IllegalArgumentException("Can't pivot in initial pivot"))
                pivot(pivotRow, j)
                canPivotRow[pivotRow] = false
            } else {
                canPivotRow[basicRowsByCol[j]] = false
            }
        }
        pivotOutAllArtificialVariables()
        columnsToInclude.forEach { assert( isBasicColumn(it) ) }
        assert(basicColsByRow.all { it != -1 })
    }


    // pivots in any rows that aren't already pivoted in
    // using a greedy algorithm:
    // for each row in turn, choose to pivot on the column with smallest
    // sparse size among all that result in a non-negative pivot value
    fun pivotOutAllArtificialVariables() {
        for(pivotRow in basicColsByRow.indices) {
            if(basicColsByRow[pivotRow] == -1) {
                val Bi = B[pivotRow]
                val iter = rows[pivotRow].iterator()
                if(iter.hasNext()) iter.advance()
                while(iter.hasNext() && Bi/iter.value() <= 0.0) {
                    iter.advance()
                }
                if(Bi/iter.value() >= 0.0) {
                    pivot(pivotRow, iter.key())
                } else {
                    throw(RuntimeException("No elements to pivot on."))
                }
            }
        }
    }



    // Does the actual pivoting of M
    // For a pivot at point (i,j), the k'th row of M and B is updated according to
    // M_k' = (M_ij*M_k - M_kj*M_i)/G
    // where G is the greatest common divisor of M_ij and M_kj
    // if the pivot point is -ve, multiplies the pivot row by -1 to make it positive
    // fun pivot(point: PivotPoint) = pivot(point.row, point.col)
    fun pivot(i: Int, j: Int) {
        var Mij = M[i,j]

        val rowIndices = rowKeys(i)
        val rowValues = Array(rowIndices.size) { m ->
            val k = rowIndices[m]
            val newRowVal = M[i,k] / Mij
            M[i,k] = newRowVal
            newRowVal
        }

        for(row in 0 until nConstraints) {
            val rowWeight = M[row,j]
            if(row != i && rowWeight != 0.0) {
                for(rowEntry in rowIndices.indices) {
                    M[row, rowIndices[rowEntry]] -= rowValues[rowEntry] * rowWeight
                }
            }
        }

        val leavingColumn = basicColsByRow[i]
        if(leavingColumn != -1) basicRowsByCol[leavingColumn] = -1
        basicColsByRow[i] = j
        basicRowsByCol[j] = i
    }


    fun rowKeys(i: Int): IntArray {
        val iter = rows[i].iterator()
        return IntArray(rows[i].size()) {
            iter.advance()
            iter.key()
        }
    }
    // Returns the rows in column j that are
    // pivot points that maintain a positive solution
    fun pivotableRows(j: Int, allowPivotsOnNegativeElements: Boolean): List<Int> {
        if(isBasicColumn(j)) return emptyList()

        var dXjmax: Double? = null
        val limits = ArrayList<Pair<Int,Double>>()
        for(i in 0 until nConstraints) {
            val Mij = M[i,j]
            if(Mij > 0.0 || (allowPivotsOnNegativeElements && B[i] <= 0.0)) {
                val dXji = B[i]/Mij
                if (dXji <= dXjmax?:dXji) dXjmax = dXji
                limits.add(Pair(i, dXji))
            }
        }

        return limits
            .filter { it.second == dXjmax }
            .map { it.first }
    }


    // true if this column is currently a basic variable
    fun isBasicColumn(j: Int): Boolean {
        return basicRowsByCol[j] != -1
    }


    interface Vector {
        operator fun get(i: Int): Double
        operator fun set(i: Int, value: Double)
    }

    interface Matrix {
        operator fun get(i: Int, j: Int): Double
        operator fun set(i: Int, j: Int, value: Double)
    }

    data class PivotPoint(val row: Int, val col: Int)

}