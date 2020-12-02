import lib.sparseIntMatrix.HashRowColIntMatrix
import lib.sparseIntMatrix.SparseIntMatrix
import lib.sparseIntMatrix.SparseIntVector
import org.apache.commons.math3.fraction.Fraction
import org.apache.commons.math3.util.ArithmeticUtils
import org.apache.commons.math3.util.Pair
import java.lang.RuntimeException
import kotlin.math.absoluteValue
import kotlin.math.sign

// Simplex algorithm with integer coefficients
class IntSimplex {
    var M: HashRowColIntMatrix            // The tableau, last row is the objective, last column is B
    val basicColsByRow: IntArray          // gives the column that is currently the basic variable for this row (-1 if none)
    val B: SparseIntVector
        get() = M.columns[bColumn]
    val objective: SparseIntVector
        get() = M.rows[objectiveRow]
    val objectiveRow: Int
        get() = M.nRows-1
    val bColumn: Int
        get() = M.nCols-1
    val X: Map<Int,Fraction>
        get() {
            val x = HashMap<Int,Fraction>(M.nRows)
            for(i in 0 until objectiveRow) {
                val basicCol = basicColsByRow[i]
                x[basicCol] = Fraction(B[i],M[i, basicCol])
            }
            return x
        }

    data class PivotPoint(val row: Int, val col: Int)

    constructor(xCoefficients: SparseIntMatrix, constants: SparseIntVector, objective: SparseIntVector) {
        val initialSolution = xCoefficients.IPsolve(constants, DoubleArray(xCoefficients.nCols) {0.0}.asList(), "==").toSparseIntVector()

        M = HashRowColIntMatrix(xCoefficients.nRows+1, xCoefficients.nCols+1)
        xCoefficients.entries.forEach { entry ->
            M[entry.row,entry.col] = entry.value
        }
        objective.forEach { (j,Cj) ->
            M[objectiveRow, j] = Cj
        }
        constants.forEach { (i, Bi) ->
            M[i,bColumn] = Bi
        }
        basicColsByRow = IntArray(xCoefficients.nRows) { -1 }
        pivotToSolution(initialSolution)
    }


    // Takes a canonical form tableaux, where the final row is the objective
    // the final column is B and the rest is M, the matrix to solve.
    // M should already have at least one column pivoted-in on each row
    constructor(canonicalFormTableaux: SparseIntMatrix) {
        M = HashRowColIntMatrix(canonicalFormTableaux)
        for((i, Bi) in B) { // ensure B is +ve
            if(Bi < 0) M.timesAssignRow(i, -1)
        }
        basicColsByRow = IntArray(M.nRows-1) { -1 }
        for(j in 0 until bColumn) {
            val col = M.columns[j]
            if(col.sparseSize == 1) {
                val (i, Mij) = col.first()
                if(Mij > 0) basicColsByRow[i] = j
            }
        }
        if(!basicColsByRow.none { it == -1 }) throw(RuntimeException("Tableaux is not in canonical form"))
    }


    // perform simplex algorithm with Bland(77) ordering
    // to find the solution that minimises the objective
    // returns true if the minimum was found. Otherwise
    // the solution is unbounded.
    fun minimise(): Boolean {
//        println("Initial state")
//        println(M)
//        var nFractional = 0
//        var nIterations = 0
        while(bland77Pivot()) {
 //           println("Pivoting...")
//            println(M)
//            if(!X.isInteger()) nFractional++
//            nIterations++
//            println("Fractional proportion = ${nFractional*1.0/nIterations}")
        }
        return objective.none { it.value < 0 }
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
        val pivotCol = objective.asSequence()
            .mapNotNull { if(it.value < 0 && it.key != bColumn) it.key else null }
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
    fun pivotToSolution(solution: SparseIntVector) {
        for((j,_) in solution) {
            if(!isBasicColumn(j)) {
                val pivotRow =
                    M.columns[j]
                        .filter { it.key != objectiveRow && basicColsByRow[it.key] == -1 }
                        .minBy { it.value.absoluteValue }
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
                val pivotCol = M.rows[pivotRow]
                    .minBy { it.value.absoluteValue*1024 + M.columns[it.key].sparseSize }
                    ?.key
                    ?:throw(RuntimeException("No elements to pivot on, equations must be degenerate"))
                pivot(pivotRow, pivotCol)
            }
        }
    }

    // true if this column is currently a basic variable
    fun isBasicColumn(j: Int): Boolean {
        val col = M.columns[j]
        if(col.sparseSize != 1) return false
        return basicColsByRow[col.first().key] == j
    }


    // Does the actual pivoting of M
    // For a pivot at point (i,j), the k'th row of M and B is updated according to
    // M_k' = (M_ij*M_k - M_kj*M_i)/G
    // where G is the greatest common divisor of M_ij and M_kj
    // if the pivot point is -ve, multiplies the pivot row by -1 to make it positive
    fun pivot(i: Int, j: Int) {
        assert(i>=0 && i < M.nRows-1)
        assert(j>=0 && j < M.nCols-1)
        var Mij = M[i,j]
        assert(Mij != 0)
        val rGCD = rowGCD(i) * Mij.sign
        if(rGCD != 1) {
            for((k, Mik) in M.rows[i]) {
                M[i,k] = Mik/rGCD
            }
            Mij /= rGCD
        }

        if(Mij.absoluteValue > 16) println("Mij is $Mij")

        val col = M.columns[j].filter { it.key != i } // separate to prevent concurrent modification
        for ((k, Mkj) in col) {
            val gcd = ArithmeticUtils.gcd(Mij, Mkj)
            M.timesAssignRow(k, Mij/gcd)
            M.weightedRowPlusAssign(k, i, -Mkj/gcd)
            rowDivideByGCD(k) // TODO: Test
        }
        basicColsByRow[i] = j
    }


    // Divides a row through by the GCD of all elements in the row
    private fun rowDivideByGCD(i: Int) {
        val rGCD = rowGCD(i)
        if(rGCD != 1) {
            for((k, Mik) in M.rows[i]) {
                M[i,k] = Mik/rGCD
            }
        }
    }

    // Calculates the GCD of all elements in the given row
    private fun rowGCD(i: Int): Int {
        var rgcd: Int = 0
        for((_, Mij) in M.rows[i]) {
            rgcd = ArithmeticUtils.gcd(rgcd, Mij)
            if(rgcd == 1) return 1
        }
        return rgcd
    }


    // Returns the rows in column j that are
    // pivot points that maintain a positive solution
    fun pivotableRows(j: Int): List<Int> {
        if(isBasicColumn(j)) return emptyList()

        var dXjmax = Fraction(Int.MAX_VALUE,1)
        val limits = M.columns[j].mapNotNull { (i, Mij) ->
            if(i != objectiveRow && Mij > 0) {
                val dXji = Fraction(B[i], Mij)
                if (dXji < dXjmax) dXjmax = dXji
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


    // True if the current solution is integer
    fun isAtIntegerSolution(): Boolean {
        return basicColsByRow.withIndex().all { entry ->
            B[entry.index].rem(M[entry.index,entry.value]) == 0
        }
    }

}