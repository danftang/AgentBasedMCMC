import lib.sparseMatrix.GridMapMatrix
import lib.sparseVector.SparseVector
import org.apache.commons.math3.linear.OpenMapRealVector
import java.lang.RuntimeException
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashMap
import kotlin.collections.HashSet
import kotlin.math.absoluteValue

class SimplifiedSimplex {
    val columns = ArrayList<TreeMap<Int, Double>>()
    val rows = ArrayList<HashSet<Int>>()
    var firstSlackColumn: Int
    val basicColsByRow: IntArray
    val epsilon = 1e-6

    val nConstraints: Int
        get() = rows.size-1

    val nVariables: Int
        get() = columns.size-1


    val M: Matrix = object: Matrix {
        override fun get(i: Int, j: Int): Double {
//            assert(i>=0)
//            assert(i<rows.size)
//            assert(j>=0)
//            assert(j<columns.size)
            return columns[j][i]?:0.0
        }

        override fun set(i: Int, j: Int, value: Double) {
            if(value.absoluteValue > epsilon) {
                if (columns[j].put(i, value) == null) rows[i].add(j)
            } else {
                if(columns[j].remove(i) != null) rows[i].remove(j)
            }
        }
    }


    val B = object: Vector {
        override fun get(i: Int) = columns.last()[i]?:0.0

        override fun set(i: Int, value: Double) {
            if(value.absoluteValue > epsilon) {
                if (columns.last().put(i, value) == null) rows[i].add(columns.lastIndex)
            } else {
                if(columns.last().remove(i) != null) rows[i].remove(columns.lastIndex)
            }
        }
    }

    val objective = object: Vector {
        override fun get(i: Int) = columns[i][rows.lastIndex]?:0.0

        override fun set(i: Int, value: Double) {
            if(value.absoluteValue > epsilon) {
                if (columns[i].put(rows.lastIndex, value) == null) rows.last().add(i)
            } else {
                if(columns[i].remove(rows.lastIndex) != null) rows.last().remove(i)
            }
        }
    }


    constructor(constraints: List<Constraint<Number>>, objective: Map<Int,Number>, initialSolution: Map<Int,Number>) {
        println("Transferring constraints to simplex tableau...")
        val nVariables = Integer.max(constraints.numVars(),
            objective.keys
                .max()
                ?.let { it + 1 }
                ?: 0
        )
        firstSlackColumn = nVariables
        val nSlackVars = constraints.count { it.relation != "==" }
        var nextSlackVar = nVariables
        val initialNonZeroColumns = ArrayList(initialSolution.keys)
        columns.ensureCapacity(nVariables + nSlackVars + 1)
        repeat(nVariables + nSlackVars + 1) { columns.add(TreeMap()) }
        rows.ensureCapacity(constraints.size + 1)
        repeat(constraints.size + 1) { rows.add(HashSet()) }
        basicColsByRow = IntArray(rows.size-1) { -1 }
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

        println("Pivoting in initial solution")
        initialPivot(initialNonZeroColumns)
        println("Done")
    }


    fun X(includeSlacks: Boolean=false): Map<Int,Double> {
        val x = HashMap<Int,Double>()
        for (i in basicColsByRow.indices) {
            val basicCol = basicColsByRow[i]
            if(includeSlacks || basicCol < firstSlackColumn) x[basicCol] = B[i] / M[i, basicCol]
        }
        return x
    }


    // pivots in all columns in the supplied list
    // then greedily pivots in any remaining rows
    fun initialPivot(columnsToInclude: List<Int>) {
        val canPivotRow = Array(basicColsByRow.size) { true }
        for(j in columnsToInclude) {
            if(!isBasicColumn(j)) {
                val pivotRow =
                    columns[j]
                        .filter { it.key != rows.lastIndex && canPivotRow[it.key] }
                        .minBy { rows[it.key].size }
                        ?.key
                        ?:throw(RuntimeException("Can't pivot to solution. Must not be an extreme solution."))
                pivot(pivotRow, j)
                canPivotRow[pivotRow] = false
            } else {
                canPivotRow[columns[j].keys.find { it != rows.lastIndex }!!] = false
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
                val pivotCol = rows[pivotRow]
                    .filter { it != columns.lastIndex &&  Bi/M[pivotRow,it] >= 0.0 }
                    .minBy { columns[it].size }
                    ?:throw(RuntimeException("No elements to pivot on."))
                pivot(pivotRow, pivotCol)
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

        val rowIndices = rows[i].toIntArray()
        val rowValues = Array(rowIndices.size) { m ->
            val k = rowIndices[m]
            val newRowVal = M[i,k] / Mij
            M[i,k] = newRowVal
            newRowVal
        }

        val colEntries = ArrayList<AbstractMap.SimpleEntry<Int,Double>>(columns[j].size - 1)
        for(entry in columns[j]) {
            if(entry.key != i) colEntries.add(AbstractMap.SimpleEntry(entry.key, entry.value))
        }

        for(rowEntry in rowValues.indices) {
            for(colEntry in colEntries) {
                M[colEntry.key, rowIndices[rowEntry]] -= rowValues[rowEntry] * colEntry.value
            }
        }
        basicColsByRow[i] = j
    }


    // Returns the rows in column j that are
    // pivot points that maintain a positive solution
    fun pivotableRows(j: Int, allowPivotsOnNegativeElements: Boolean): List<Int> {
        if(isBasicColumn(j)) return emptyList()

        var dXjmax: Double? = null
        val limits = columns[j].mapNotNull { (i, Mij) ->
            if(i != rows.lastIndex && (Mij > 0.0 || (allowPivotsOnNegativeElements && B[i] <= 0.0))) {
                val dXji = B[i]/Mij
                if (dXji <= dXjmax?:dXji) dXjmax = dXji
                Pair(i, dXji)
            } else null
        }

        return limits
            .filter { it.second == dXjmax }
            .map { it.first }
    }


    // true if this column is currently a basic variable
    fun isBasicColumn(j: Int): Boolean {
        val col = columns[j]
        return if (col.size > 2) {
            false
        } else {
            col.keys
                .find { it != rows.lastIndex }
                ?.let { basicColsByRow[it] == j }
                ?:false
        }
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