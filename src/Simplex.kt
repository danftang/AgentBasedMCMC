import com.google.ortools.linearsolver.MPSolver
import lib.abstractAlgebra.*
import lib.sparseMatrix.GridMapMatrix
import lib.vector.MapVector
import lib.vector.MutableMapVector
import java.lang.RuntimeException


class Simplex<T: Any>(
    val M: GridMapMatrix<T>,
    val basicColsByRow: IntArray = IntArray(M.nRows-1) { -1 }
): FieldOperators<T> by M {

    val B: MutableMapVector<T>
        inline get() = M.columns[bColumn]
    val objective: MutableMapVector<T>
        inline get() = M.rows[objectiveRow]
    val objectiveRow: Int
        inline get() = M.nRows-1
    val bColumn: Int
        inline get() = M.nCols-1

    data class PivotPoint(val row: Int, val col: Int)


    fun X(): Map<Int,T> {
        val x = HashMap<Int, T>()
        for (i in 0 until objectiveRow) {
            val basicCol = basicColsByRow[i]
            x[basicCol] = B[i] / M[i, basicCol]!!
        }
        return x
    }


    // perform simplex algorithm with Bland(77) ordering
    // to find the solution that minimises the objective
    // returns true if the minimum was found. Otherwise
    // the solution is unbounded.
    fun minimise(): Boolean {
        while (bland77Pivot()) {
            println(M)
        }
        return objective.nonZeroEntries.none { it.value < zero }
    }


    // Perform the next pivot according to the ordering
    // defined in the first algorithm in Bland(77)
    // i.e. choose the column with the lowest index from among those
    // that have -ve coefficient in the objective.
    // If no columns have -ve coefficient, we have reached a minimum.
    // Then choose the row whose basic variable has the lowest (column) index
    // from among the columns with +ve coefficient.
    // If no columns have +ve coefficient, then the solution is unbounded
    // Returns true if the pivot was performed
    private fun bland77Pivot(): Boolean {
        val pivotCol = objective.nonZeroEntries.asSequence()
            .mapNotNull { if(it.value < zero && it.key != bColumn) it.key else null }
            .min() // TODO: Could be made faster by keeping the objective row in an ordered tree
            ?:return false
        val pivotRow = pivotableRows(pivotCol)
            .minBy { basicColsByRow[it] }
            ?:return false
        pivot(pivotRow, pivotCol)
        return true
    }


    // pivots in all non-zero columns in the supplied solution
    // then greedily pivots in any remaining rows
    fun pivotToSolution(solution: Map<Int,T>) {
        for((j,_) in solution) {
            if(!isBasicColumn(j)) {
                val pivotRow =
                    M.columns[j].nonZeroEntries
                        .filter { it.key != objectiveRow && basicColsByRow[it.key] == -1 }
                        .minBy { M.rows[it.key].nonZeroEntries.size }
//                        .minBy { it.value.absoluteValue }
                        ?.key
                        ?:throw(RuntimeException("Can't pivot to solution. Must be a non-extreme solution (i.e. contains loops)!"))
                pivot(pivotRow, j)
            }
        }
        greedyPivotAll()
    }


    // pivots in any rows that aren't already pivoted in
    // using a greedy algorithm:
    // for each row in turn, choose to pivot on the column with smallest
    // sparse size among all that have a +-1 in this row
    fun greedyPivotAll() {
        for(pivotRow in basicColsByRow.indices) {
            if(basicColsByRow[pivotRow] == -1) {
                val pivotCol = M.rows[pivotRow].nonZeroEntries
                    .minBy { M.columns[it.key].nonZeroEntries.size }
//                    .minBy { it.value.absoluteValue*1024 + M.columns[it.key].size }
                    ?.key
                    ?:throw(RuntimeException("No elements to pivot on, equations must be degenerate"))
                pivot(pivotRow, pivotCol)
            }
        }
    }

    // true if this column is currently a basic variable
    fun isBasicColumn(j: Int): Boolean {
        val col = M.columns[j]
        if(col.nonZeroEntries.size != 1) return false
        return basicColsByRow[col.nonZeroEntries.entries.first().key] == j
    }


    // Does the actual pivoting of M
    // For a pivot at point (i,j), the k'th row of M and B is updated according to
    // M_k' = (M_ij*M_k - M_kj*M_i)/G
    // where G is the greatest common divisor of M_ij and M_kj
    // if the pivot point is -ve, multiplies the pivot row by -1 to make it positive
    fun pivot(i: Int, j: Int) {
        assert(i>=0 && i < M.nRows-1)
        assert(j>=0 && j < M.nCols-1)
        var Mij = M[i,j]!!

        M.rows[i] /= Mij
//        M.rowReassign(i) { it/Mij }
        val colEntries = M.columns[j].nonZeroEntries.asSequence().filter { it.key != i }.toList()
        val rowEntries = M.rows[i].nonZeroEntries.entries.toList()
        for(rowEntry in rowEntries) {
            for(colEntry in colEntries) {
                val outerProd = rowEntry.value * colEntry.value
                M.mapAssign(colEntry.key, rowEntry.key) { oldVal -> oldVal - outerProd }
//                M.compute(colEntry.key, rowEntry.key) { _, Mpq ->
//                    val newVal = (Mpq?:zero) - outerProd
//                    if(newVal == zero) null else newVal
//                }
            }
        }
        basicColsByRow[i] = j
    }

    // Returns the rows in column j that are
    // pivot points that maintain a positive solution
    fun pivotableRows(j: Int): List<Int> {
        if(isBasicColumn(j)) return emptyList()

        var dXjmax: T? = null
        val limits = M.columns[j].nonZeroEntries.mapNotNull { (i, Mij) ->
            if(i != objectiveRow && Mij > zero) {
                val dXji = B[i]/Mij
                if (dXji <= dXjmax?:dXji) dXjmax = dXji
                Pair(i, dXji)
            } else null
        }

        return limits
            .filter { it.second == dXjmax }
            .map { it.first }
    }


    fun allPositivePivotPoints(): List<PivotPoint> {
        return (0 until bColumn).asSequence()
            .flatMap { j ->
                pivotableRows(j).asSequence()
                    .map { i ->
                        PivotPoint(i,j)
                    }
            }.toList()
    }

}


fun MPSolver.toSimplex(): Simplex<Double> {
    val grid = GridMapMatrix<Double>(DoubleOperators,numConstraints()+1, numVariables()+1)

    // check format
    constraints().forEach { constraint ->
        if(constraint.lb() != constraint.ub()) throw(RuntimeException("Can only do equality constraints at the moment"))
    }
    variables().forEach { variable ->
        if(variable.lb() != 0.0 || variable.ub() != Double.POSITIVE_INFINITY) throw(RuntimeException("Can only do positive variables for now"))
    }

    constraints().forEach { constraint ->
        variables().forEach { variable ->
            val coeff = constraint.getCoefficient(variable)
            if(coeff != 0.0) grid[constraint.index(), variable.index()] = coeff
        }
        grid[constraint.index(), grid.nCols-1] = constraint.lb()
    }

    variables().forEach { variable ->
        val coeff = objective().getCoefficient(variable)
        grid[grid.nRows-1, variable.index()] = coeff
    }

    val simplex = Simplex(grid)

    this.solve()
    val initialSolution = variables()
        .filter { it.solutionValue() != 0.0 }
        .associate { Pair(it.index(), it.solutionValue()) }

    simplex.pivotToSolution(initialSolution)

    return simplex
}
