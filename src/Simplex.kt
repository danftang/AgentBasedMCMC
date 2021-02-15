import com.google.ortools.linearsolver.MPSolver
import lib.abstractAlgebra.*
import lib.sparseMatrix.GridMapMatrix
import lib.sparseMatrix.IPsolve
import lib.sparseMatrix.SparseMatrix
import lib.sparseMatrix.copyTo
import lib.vector.*
import org.apache.commons.math3.fraction.Fraction
import java.lang.Integer.max
import java.lang.RuntimeException
import java.util.AbstractMap
import kotlin.math.roundToInt

//////////////////////////////////////////////////////////////////////////
// Represents a linear program in the form
// maximise CX
// subject to
// AX = B
// X >= 0
//
// Packed into a single matrix
//
// M = (A|B)
//     (C 0)
//////////////////////////////////////////////////////////////////////////
open class Simplex<T>(
    val M: GridMapMatrix<T>,
    val basicColsByRow: IntArray = IntArray(M.nRows-1) { -1 }
): FieldOperators<T> by M
    where T: Number,
          T: Comparable<T>
{

    val B: MutableSparseVector<T>
        inline get() = M.columns[bColumn]
    val objective: MutableSparseVector<T>
        inline get() = M.rows[objectiveRow]
    val objectiveRow: Int
        inline get() = M.nRows-1
    val bColumn: Int
        inline get() = M.nCols-1

    data class PivotPoint(val row: Int, val col: Int)


    fun X(): SparseVector<T> {
        val x = B.new()
        for (i in 0 until objectiveRow) {
            val basicCol = basicColsByRow[i]
            x[basicCol] = B[i] / M[i, basicCol]
        }
        return x
    }

    constructor(
        xCoefficients: SparseMatrix<T>,
        constants: SparseVector<T>,
        objective: SparseVector<T>)
            : this(
        GridMapMatrix(xCoefficients.operators,xCoefficients.nRows+1, xCoefficients.nCols+1)
    ) {
//        val initialSolution = xCoefficients.IPsolve(constants, DoubleArray(xCoefficients.nCols) {0.0}.asList(), "==")
        xCoefficients.copyTo(M)
        objective.copyTo(M.rows[objectiveRow])
        constants.copyTo(M.columns[bColumn])

        // Find initial positive solution
        val initialSolution = xCoefficients
            .IPsolve(constants, emptyMap<Int,T>().asVector(xCoefficients.operators), "==")
            .asList()
            .mapIndexedNotNull { index, value ->
                if(value != 0.0) index else null
            }

        pivotColumns(initialSolution)
    }

    constructor(constraints: List<Constraint<T>>, objective: SparseVector<T>) : this(
        GridMapMatrix(objective.operators, constraints.size+1, 1)
    ) {
        val nVariables = max(
            constraints
                .flatMap { it.coefficients.keys }
                .max()
                ?:-1,
            objective.nonZeroEntries.keys
                .max()
                ?:-1
        ) + 1
        val nSlackVars = constraints.count { it.relation != "==" }
        var nextSlackVar = nVariables
        M.resize(M.nRows, nVariables + nSlackVars + 1)
        constraints.forEachIndexed { i, constraint ->
            constraint.coefficients.forEach { (j, x) ->
                M[i,j] = x
            }
            when(constraint.relation) {
                ">=" -> M[i,nextSlackVar++] = -operators.one
                "<=" -> M[i,nextSlackVar++] = operators.one
            }
            B[i] = constraint.constant
        }
        objective.nonZeroEntries.forEach { (j,x) ->
            this.objective[j] = x
        }

        println("Constraints in M")
        println(M)
        println("Finding initial solution")
        findInitialSolution()

    }


    // Find an initial positive pivot state from a raw set of constraints.
    //
    // First look for any columns that are already potentially basic
    // then pivot in the remaining rows on the first element
    // then pivot to reduce negativity until non-negative solution is found
    fun findInitialSolution() {
        for(j in 0 until bColumn) {
            if(M.columns[j].nonZeroEntries.size == 1) {
                val i = M.columns[j].nonZeroEntries.keys.first()
                basicColsByRow[i] = j
                M.rows[i] /= M[i,j]
            } else if (M.columns[j].nonZeroEntries.size == 2 && !objective[j].isZero()) {
                val i = M.columns[j].nonZeroEntries.keys.min()!!
                basicColsByRow[i] = j
                M.rows[i] /= M[i,j]
            }
        }
        for(i in basicColsByRow.indices) {
            if(basicColsByRow[i] < 0) {
               pivot(i, M.rows[i].nonZeroEntries.keys.first())
            }
        }
        println("Initial pivoted solution")
        println(M)
        pivotOutNegatives()
    }

    // Pivots until all elements of X are non-negative.
    //
    // Does this by repeatedly pivoting on the column
    // whose sum over elements in rows with negative B is minimised
    // until either X is non-negative or no column has a negative sum
    // over elements in B-negative rows.
    // This is guaranteed to get to a positive solution or prove that
    // none exists since if no column has a negative sum then no pivot increases the
    // sum of negative elements of X. So the sum of these elements of X must be maximised
    // (by the Simplex algorithm), but if this is true at a negative solution then
    // there must not exist a non-negative solution since this would have a higher
    // sum than the maximum.
    //
    // returns true if a non-negative solution is found.
    fun pivotOutNegatives(): Boolean {
        println("Pivoting out negatives")
        var negativeRows: List<MutableMapVector<T>>
        do {
            negativeRows = B.nonZeroEntries
                .entries
                .filter { it.value < operators.zero && it.key < objectiveRow }
                .map { M.rows[it.key] }
            val sumOfNegativeRows = MutableMapVector(operators, HashMap())
            negativeRows.forEach {
                sumOfNegativeRows += it
            }
            sumOfNegativeRows.nonZeroEntries.remove(bColumn)
            val (pivotColumn, bNegativeSum) = sumOfNegativeRows.nonZeroEntries
                .minBy { it.value }
                ?: AbstractMap.SimpleEntry<Int, T>(0, operators.zero)
            if (bNegativeSum < operators.zero) {
                // find row to pivot on
                var pivotRow: Pair<Int,T>? = null
                M.columns[pivotColumn].nonZeroEntries.forEach { (i, Mij) ->
                    val dXj = B[i]/Mij
                    if(i != objectiveRow && dXj > zero && (pivotRow?.let { dXj < it.second } != false)) {
                        pivotRow = Pair(i, dXj)
                    }
                }
                pivot(pivotRow!!.first, pivotColumn) // pivotRow cannot be null if bNegativeSum < 0
                println(M)
                println()
            }
        } while(bNegativeSum < operators.zero)
        return negativeRows.isEmpty()
    }

    // perform simplex algorithm with Bland(77) ordering
    // to find the solution that minimises the objective
    // returns true if the minimum was found. Otherwise
    // the solution is unbounded.
    fun minimise(): Boolean {
        while (bland77Pivot()) {
            println(M)
        }
        return objective.nonZeroEntries.none { it.value < -zero }
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
            .mapNotNull { if(it.value < -zero && it.key != bColumn) it.key else null }
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
    fun pivotColumns(columnsToPivotIn: List<Int>) {
        for(j in columnsToPivotIn) {
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
        val row = col.nonZeroEntries.entries.first().key
        if(row == objectiveRow) return false
        return basicColsByRow[row] == j
    }


    // Does the actual pivoting of M
    // For a pivot at point (i,j), the k'th row of M and B is updated according to
    // M_k' = (M_ij*M_k - M_kj*M_i)/G
    // where G is the greatest common divisor of M_ij and M_kj
    // if the pivot point is -ve, multiplies the pivot row by -1 to make it positive
    inline fun pivot(point: PivotPoint) = pivot(point.row, point.col)
    fun pivot(i: Int, j: Int) {
        assert(i>=0 && i < M.nRows-1)
        assert(j>=0 && j < M.nCols-1)
        var Mij = M[i,j]

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
            .filter { it.second as Number == dXjmax }
            .map { it.first }
    }


    // maximum +ve value this column can take in a pivot without forcing -ve solution
    // null if there is no limit.
    fun columnPivotLimit(j: Int): T? {
        var dXjmax: T? = null
        M.columns[j].nonZeroEntries.forEach { (i, Mij) ->
            if(i != objectiveRow && Mij > zero) {
                val dXji = B[i]/Mij
                if (dXji <= dXjmax?:dXji) dXjmax = dXji
            }
        }
        return dXjmax
    }


    fun allPositivePivotPoints(): List<PivotPoint> {
        return (0 until M.nCols).asSequence()
            .filter { it != bColumn }
            .flatMap { j ->
                pivotableRows(j).asSequence()
                    .map { i ->
                        PivotPoint(i,j)
                    }
            }.toList()
    }

//    data class PivotPointsByDegeneracy(val degeneratePivots: List<PivotPoint>, val nonDegeneratePivots: List<PivotPoint>)
//    fun allPositivePivotPointsByDegeneracy(): PivotPointsByDegeneracy {
//        val degeneratePivots = ArrayList<PivotPoint>()
//        val nonDegeneratePivots = ArrayList<PivotPoint>()
////        return (0 until M.nCols) asSequence()
////            .filter { it != bColumn }
////            .flatMap { j ->
////                pivotableRows(j).asSequence()
////                    .map { i ->
////                        PivotPoint(i,j)
////                    }
////            }.toList()
//    }

    fun isDegenerate(pivot: PivotPoint): Boolean = B[pivot.row].isZero()


    // returns all pivots that maintain a positive solution
    // and that do not have a zero in the B column
    // (i.e. that a pivot here will change the solution)
    fun allNonDegeneratePivots(): List<PivotPoint> {
        return (0 until M.nCols).asSequence()
            .filter { it != bColumn }
            .flatMap { j ->
                pivotableRows(j).asSequence()
                    .filter { i -> !B[i].isZero() }
                    .map { i ->
                        PivotPoint(i,j)
                    }
            }.toList()

    }

    // returns all pivots that maintain a positive solution
    // and that have a zero in the B column
    // (i.e. that a pivot here will not change the solution)
    fun allDegeneratePivots(): List<PivotPoint> {
        return (0 until bColumn).asSequence()
            .flatMap { j ->
                pivotableRows(j).asSequence()
                    .filter { i -> B[i].isZero() }
                    .map { i ->
                        PivotPoint(i,j)
                    }
            }.toList()
    }




    fun postPivotB(pivot: PivotPoint): SparseVector<T> {
        val Bpiv = B[pivot.row]/M[pivot.row,pivot.col]
        val Bprime = B.toMutableSparseVector()
        Bprime.weightedPlusAssign(M.columns[pivot.col], -Bpiv)
        Bprime[pivot.row] = Bpiv
        return Bprime
    }



    fun T.roundToInt(): Int {
        return when(this) {
            is Double -> roundToInt()
            is Fraction -> (numerator*1.0/denominator).roundToInt()
            else -> throw(UnsupportedOperationException("Don't know how to round a ${this.javaClass.simpleName}"))
        }
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
        .mapNotNull { if(it.solutionValue() != 0.0) it.index() else null }

    simplex.pivotColumns(initialSolution)

    return simplex
}
