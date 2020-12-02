import lib.sparseIntMatrix.HashIntVector
import lib.sparseIntMatrix.HashRowColIntMatrix
import lib.sparseIntMatrix.IntArrayVector
import lib.sparseIntMatrix.SparseIntVector
import kotlin.math.sign

// Represents the convex polyhedron defined by the points X such that
// MX = B
// and for all elements x_i of X:
// x_i >= 0
// x_i is integer
//
// We allow various linear transformations of the polyhedron so that, more
// generally
// VM(UX' + R) = VB
// such that MX = B <-> VMU(X' + R) = VB where X = U(X' + R)
// and both X and X' are in the non-negative integer space.
//
// If TMX' = TB
class ConvexPolyhedron(var M: HashRowColIntMatrix, var B: HashIntVector) {
    var U = HashRowColIntMatrix.identity(M.nCols) // Linear transform of
    var R = HashIntVector()

    // find the valid point that maximises objective.transpose()*X
    // throws RuntimeException if there are no valid points in this polyhedron
    fun maximise(objective: List<Double>): IntArrayVector {
        return U * M.IPsolve(B, objective, "=") + R
    }

    // returns a random valid point
    fun findValidPoint(): IntArrayVector {
        return maximise(DoubleArray(M.nCols) {0.0}.asList())
    }


    // removes solution dimensions which can only be zero
    // i.e. all unknowns that appear in rows whose coefficients are all one sign
    fun removeZeros() {
        // first find redundant rows and cols
        // val T = SparseColIntMatrix.Identity(M.nCols)
        var isReduced: Boolean
        do {
            isReduced = false
            val redundantCols = HashSet<Int>() // cols with appearing in rows with all same-sign coefficients
            val redundantRows = HashSet<Int>() // rows with all zero coefficients
            for(i in 0 until M.nRows) {
                if(B[i] == 0) {
                    val row = M.rows[i]
                    val rowType = row.values.fold<Int, Int?>(0) { acc, v ->
                        if (acc == null || acc * v.sign == -1) null else v.sign
                    }
                    when (rowType) {
                        -1, 1 -> redundantCols.addAll(row.keys)
                        0 -> redundantRows.add(i)
                    }
                }
            }

            // now remove redundant cols and rows
            if(!redundantRows.isEmpty()) {
                val R = HashRowColIntMatrix.identity(M.nRows)
                R.removeRows(redundantRows)
                M = R * M
                B = R * B
                isReduced = true
            }
            if(!redundantCols.isEmpty()) {
                U.removeColumns(redundantCols)
                M.removeColumns(redundantCols)
                isReduced = true
            }
        }while(isReduced)
    }

    // constrains the solution so that certain
    // elements of the solution have certain values
    // removing the columns that correspond to the constrained values
    // and updating B and R appropriately
    // constraints specifies a map from solution index to value for that element
    fun constrainSolution(constraints: SparseIntVector) {
        R.plusAssign(U*constraints)
        B.minusAssign(M*constraints)
        M.removeColumns(constraints.keys)
        U.removeColumns(constraints.keys)
    }
}