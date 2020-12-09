package experiments

import ABMMatrices.twoDabmMatrix
import lib.isInteger
import lib.sparseIntMatrix.HashIntVector
import lib.sparseIntMatrix.IntArrayVector
import lib.unaryMinus
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import java.lang.RuntimeException
import kotlin.random.Random


// Experiments to show that we can/can't pivot between two solutions
// via all integer, positive solutions
class PivotBetweenExpts {
    init {
        System.loadLibrary("jniortools")
    }


    // Show that we can pivot between any two solutions
    // while maintaining positive, integer intermediate
    // solutions
//    @Test
//    fun pivotBetweenSolutions() {
//        // setup the problem
//        val gridSize = 8
//        val timesteps = 6
//        val nAgents = 20
//        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
//        val observations = HashIntVector()
//        val startPos = Array(nAgents) {
//            Pair(Random.nextInt(gridSize-1), Random.nextInt(gridSize-1))
//        }
//        for(agent in startPos.indices) {
//            val xPos = startPos[agent].first
//            val yPos = startPos[agent].second
//            observations[xPos*gridSize + yPos] = -1
//            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
//        }
//
//        // first find two random solutions
//        val Xa = HashIntVector(
//            abmMatrix.IPsolve(observations, DoubleArray(abmMatrix.nCols) { Random.nextDouble() }.asList(), "==")
//        )
//        val Xb = HashIntVector(
//            abmMatrix.IPsolve(observations, DoubleArray(abmMatrix.nCols) { Random.nextDouble() }.asList(), "==")
//        )
//        println(Xa)
//        println(Xb)
//        val simplex = OldSimplex(abmMatrix, initialSolution = Xa)
//        println("Pivoted to ${simplex.basicColsByRow.asList()}")
//
//        // now make path to other solution
//        while(simplex.pivotTowardsObjective(Xb)) {
//            assert(simplex.X.isPositive())
//            println("not Xb score = ${simplex.notXbScore(Xb)}")
//            println("Fraction score = ${simplex.X.fractionScore()}")
//            println("Min fraction score on one pivot = ${simplex.minPivotFractionScore()}")
////            println("")
////            println("Done pivot")
//        }
//        assert(simplex.X.roundToIntVector() == Xb)
//
//        println("Done!!!")
//    }

//    @Test
//    fun multipleTest() {
//        for(i in 1..100) pivotBetweenSolutions()
//    }


//    fun OldSimplex.pivotTowardsObjective(Xb: HashIntVector): Boolean {
//        if(X == HashFractionVector(Xb)) return false
//
////        val xScore = notXbScore(Xb)
//        val pivotCol = (0 until M.nCols)
//            .filter { j -> !isBasicColumn(j) }
//            .map { j ->
//                applyNullVector(j)
//                val newScore = notXbScore(Xb)
//                applyNullVector(j, -X[j])
////                println("Found new score $newScore")
//                Pair(j,newScore)
//            }
//            .minBy { it.second }
//            ?.first
//            ?:throw(RuntimeException("No columns reduce the objective"))
//
//        val pivotRow = pivotableRows(pivotCol).random()
//        pivotAt(OldSimplex.PivotPoint(pivotRow,pivotCol))
//        return true
//    }
//
//
//    fun OldSimplex.notXbScore(Xb: HashIntVector): Double {
//        return X.sumByDouble { if(Xb[it.key] == 0) it.value.toDouble() else 0.0 }
//    }
//
//
//    fun OldSimplex.pivotTowardsFractional(Xb: HashIntVector): Boolean {
//        val colsToBePivotedIn = Xb.keys.toMutableSet()
//        colsToBePivotedIn.removeAll(basicColsByRow.asList())
//        if(colsToBePivotedIn.isEmpty()) { // Xb already fully pivoted in
//            return false
//        }
//        println("${colsToBePivotedIn.size} cols to be pivoted in")
//
//        val pivotCol = colsToBePivotedIn
//            .find { j ->
//                pivotableRows(j).filter { i -> Xb[basicColsByRow[i]] == 0}.isNotEmpty()
//            }
//            ?:colsToBePivotedIn
//                .find { j ->
//                    pivotableRows(j).isNotEmpty()
//                }
//            ?:
//            throw(RuntimeException("Can't pivot any rows in!"))
//        val pivotRow = pivotableRows(pivotCol)
//            .filter { i -> Xb[basicColsByRow[i]] == 0}
//            .firstOrNull()
//            ?:pivotableRows(pivotCol).random()
//        pivotAt(OldSimplex.PivotPoint(pivotRow,pivotCol))
//        return true
//    }
//
//
//    fun OldSimplex.minPivotFractionScore(): Int {
//        return (0 until M.nCols)
//            .mapNotNull {  j ->
//                val perturbation = perturbationRange(j).second
//                var score: Int? = null
//                if(perturbation != Fraction.ZERO) {
//                    applyNullVector(j, perturbation)
//                    score = X.fractionScore()
//                    applyNullVector(j, -perturbation)
//                }
//                score
//            }
//            .min()
//            ?:X.fractionScore()
//    }
//

//    fun HashFractionVector.fractionScore(): Int {
//        return sumBy { if(it.value.isInteger()) 0 else 1 }
//    }

    // makes a pivot that reduces the distance to solution Xb
    // Where the distance is measured as the number of non-zero events
    // in Xb that are not pivoted in. The pivot must maintain the positive,
    // integer nature of the current solution.
    //
    // Algorithm is:
    // Let C be the set of all non-zero events in Xb that aren't currently pivoted in
    // If there is a +-1 pivot point on a column in C and a row with b=0 that won't pivot out
    // a non-zero member of Xb then pivot there. Since this doesn't change the solution, it is guaranteed
    // to remain positive and integer.
    // Otherwise, there must be at least one column in C such that each non-zero value
    // is either on a b!=0 row or on a row with pivoted value in Xb. Pivot this in on any
    // +-1 pivot point on a b!=0 row. This should always give a positive, integer solution.
//    fun Simplex.pivotTowards(Xb: HashIntVector): Boolean {
//        val pivotedIn = pivotPoints.toSet()
//        val colsToBePivotedIn = Xb.keys.filter { !pivotedIn.contains(it) }
//        if(colsToBePivotedIn.isEmpty()) return false  // Xb already fully pivoted in
//        // look for zero pivots
//        colsToBePivotedIn.forEach { j ->
//            M.columns[j]
//                .find { (i, Mij) ->
//                    Mij.absoluteValue == 1 && B[i] == 0 && Xb[pivotPoints[i]] == 0
//                }
//                ?.also { (i, _) ->
//                    println("pivoting on zero")
//                    pivotAt(Simplex.PivotPoint(i, j))
//                    return@pivotTowards true
//                }
//        }
//        // look for loop pivots
//        colsToBePivotedIn.forEach { j ->
//            if(M.columns[j].all { (i, _) -> B[i] != 0 || Xb[pivotPoints[i]] != 0 }) {
//                println("Found loop ${M.columns[j]} Xa members ${M.columns[j].filter { B[it.key] > 0 }} B = ${B.filter {M.columns[j].keys.contains(it.key)} }")
//                M.columns[j]
//                    .find { (i, Mij) -> Xb[pivotPoints[i]] == 0 && isValidPivotPoint(Simplex.PivotPoint(i, j)) }
//                    ?.also { (i, _) ->
//                        println("pivoting on loop")
//                        pivotAt(Simplex.PivotPoint(i, j))
//                        return@pivotTowards true
//                    }
////                    ?:throw(RuntimeException("Found loop with no valid pivot point ${M.columns[j]}"))
//            }
//        }
//        println("Cols to be pivoted in ${colsToBePivotedIn}")
//        throw(RuntimeException("Can't find a zero pivot or a loop pivot. This shouldn't happen"))
//    }


    // makes a pivot that reduces the distance to Xb
    // if more than one possible pivot exists
    // then choose one at random.
    // returns true if a pivot is made, false otherwise
//    fun Simplex.pivotTowardsRandom(Xb: HashIntVector): Boolean {
//        val colsToBePivotedIn = Xb.keys.toMutableSet()
//        colsToBePivotedIn.removeAll(basicColsByRow.asList())
//        if(colsToBePivotedIn.isEmpty()) { // Xb already fully pivoted in
//            return false
//        }
//        println("${colsToBePivotedIn.size} cols to be pivoted in")
//
//        if(colsToBePivotedIn.size == 2) {
//            colsToBePivotedIn.forEach { j ->
//                println()
//                val XaRows = M.columns[j].keys.filter { i -> X[basicColsByRow[i]] != Fraction.ZERO }.toSet()
//                val XbRows = M.columns[j].keys.filter { i -> Xb[basicColsByRow[i]] != 0 }.toSet()
//
//                println("Xa and Xb = ${ XaRows.intersect(XbRows).associateWith { i -> "${M[i,j]}:${X[basicColsByRow[i]]}" }.entries.sortedBy { it.key } }")
//                println("Xa not Xb = ${ XaRows.minus(XbRows).associateWith { i -> "${M[i,j]}:${X[basicColsByRow[i]]}" }.entries.sortedBy { it.key } }")
//                println("Xb not Xa = ${ XbRows.minus(XaRows).associateWith { i -> M[i,j] }.entries.sortedBy { it.key } }")
//                println("Not Xa not Xb = ${
//                    M.columns[j].filter { (i, Mij) -> !XaRows.contains(i) && !XbRows.contains(i) }.sortedBy { it.key }
//                }")
//            }
//            val intersection = M.columns[colsToBePivotedIn.first()].keys.intersect(M.columns[colsToBePivotedIn.last()].keys)
//            println("Loop intersection ${
//                intersection
//                    .associateWith { Triple(M[it, colsToBePivotedIn.first()], M[it, colsToBePivotedIn.last()], X[basicColsByRow[it]]) }
//                    .entries
//                    .sortedBy { it.key }
//            }")
//        }
//
//        var validPivots = colsToBePivotedIn
//            .flatMap { j ->
//                M.columns[j].keys
//                    .map { i -> Simplex.PivotPoint(i,j) }
//                    .filter { !Xb.keys.contains(pivotPoints[it.row]) && isValidPivotPoint(it) }
//            }
//            .ifEmpty {
//                // need to pivot out a member of Xb
//                println("Choosing random pivot *****************************")
////                colsToBePivotedIn
////                    .flatMap {  j ->
////                        M.columns[j].keys
////                            .map { i -> Simplex.PivotPoint(i,j) }
////                            .filter { isValidPivotPoint(it) }
////                    }
//
////                findAllValidPivotPoints()
//
//                // find all pivots that reduce the pivot distance of the first col to be pivoted in
////                val col = colsToBePivotedIn[0]
////                findAllValidPivotPoints().filter { pivot -> pivotScore(pivot, col) < 0 }
//
////                findAllValidPivotPoints().filter { pivot ->
////                    colsToBePivotedIn.sumBy { col -> pivotScore(pivot, col) } < 0
////                }
//
////                findAllValidPivotPoints().filter { pivot ->
////                    colsToBePivotedIn.find { col -> pivotScore(pivot, col) < 0 } != null // pivot reduces at least one col distance
////                }
//
//                // find a pivot that zero's an element of B on a row with a unit pivot point in the cols to be pivoted-in
//                findAllValidPivotPoints().filter { pivot ->
//                    val k = B[pivot.row]/M[pivot.row,pivot.col]
//                    M.columns[pivot.col].find { (row, Mij) ->
//                        B[row] - k*Mij == 0 && colsToBePivotedIn.find { col -> M[row,col].absoluteValue == 1 } != null
//                    } != null
//                }
//
//            }
//            .ifEmpty {
//                throw(RuntimeException("No valid pivot points found."))
//            }
//        println("${validPivots.size} possible pivots")
//        val chosenPivot = validPivots.random()
//        println("Pivoting on M(${chosenPivot.row}, ${chosenPivot.col}) = ${M[chosenPivot.row, chosenPivot.col]}")
//        println("Pivoting in ${if(!Xb.keys.contains(chosenPivot.col)) "not" else ""} Xb")
//        println("Pivoting out ${if(B[chosenPivot.row] == 0) "not" else ""} Xa ${if(!Xb.keys.contains(pivotPoints[chosenPivot.row])) "not" else ""} Xb")
//        pivotAt(chosenPivot)
//        return true
//    }


//    fun Simplex.pivotAtRandomValid(): Boolean {
////        val validPivots = (0 until M.nCols).flatMap { j ->
////            M.columns[j].keys
////                .map {i -> Simplex.PivotPoint(i,j)}
////                .filter { isValidPivotPoint(it) }
////        }
//        val validPivots = findAllValidPivotPoints()
//        if(validPivots.isEmpty()) return false
//        val chosenPoint = validPivots.random()
//        println("Pivoting at $chosenPoint")
//        pivotAt(chosenPoint)
//        return true
//    }


    // Returns the change in score of a given column by pivoting at "pivot".
    // The score of a column, C, is defined as H(C-B)(C-B)
    // where H is the Heaviside step function.
    // Assumes that "pivot" is a valid pivot point
    //
    // If "col" is not the pivot column, then the change
    // comes through changes in B and C, and any sign change on the pivot row
    // If "col" is the pivot column, then we assume the pivot is
    // valid, so the current column has z pivot distance of zero, and
    // the pivoted-in column is necessarily pivotable (by the reverse pivot)
    // so the score is zero
    //
    // Let j be the column we are assessing
    // Let Mij' and Bi' be the elements of M and B after the pivot
    // Let Mpp be the value at the pivot point, Mpj be the pivot-column elements
    // and Mip be the pivot-row elements
    //
    // Mij' = Mij - Mpj*Mip/Mpp
    // Bi' = Bi - Bp*Mip/Mpp
    // So Bi' - Mij' = Bi - Bp*Mip/Mpp - Mij + Mpj*Mip/Mpp
    // = Bi - Mij + Mip(Mpj - Bp)/Mpp
    //
    // ----
    // We want to consider the pivoting-out of a member of Xb and pivoting-in of another member of Xb
    // or pivoting out and in non members of Xb that results in currently pivoted-out columns becoming
    // more pivotable. A column is pivotable if:
    //  - there exists a unit element on a row outside of Xb and Xa (i.e. Bi == 0 and not pivoting out a row in Xb)
    //  - there exists a unit element on a row in Xa and not Xb and B does not become negative after pivot
    // * Note that if Bi != 0 then Bi > 0 and we can only pivot on +ve elements
    // * If it's semi-pivotable and there exists a unit pivot not in Xb with B <=1 then it is pivotable,
    // * If Xb is an extreme point, then semi-pivotability implies pivotability (if there exists a positive, unit pivot point), so distance to semi-pivotability is a sensible measure(?)
    // * If we can always reduce the sum of semi-pivotability distances of all unpivoted columns then there exists a path
    // * If we can always reduce the semi-pivotability distance of a single column then there exists a path (this is likely easier to prove possible)
    //
//    fun Simplex.pivotScore(pivot: Simplex.PivotPoint, col: Int): Int {
//        if(col == pivot.col) return 0
//        val Mpp = M[pivot.row,pivot.col]
//        val k = (M[pivot.row,col] - B[pivot.row])/Mpp
//        return M.columns[pivot.col]
//            .map { (row, Mip) ->
//                val oldMeasure =  B[row] - M[row,col]
//                if(row != pivot.row) {
//                    negativeMagnitude(oldMeasure + Mip * k) - negativeMagnitude(oldMeasure)
//                } else {
//                    negativeMagnitude(oldMeasure/Mpp) - negativeMagnitude(oldMeasure)
//                }
//            }
//            .sum()
//    }


    // Define score for a column as the smallest negativeMagnitude possible via a unit pivot
//    fun Simplex.columnScore(col: Int): Int {
//        return M.columns[col]
//            .filter { it.value == 1 || (it.value == -1 && B[it.key] == 0) }
//            .map { (pivotRow, pivotVal) ->
//                val k = B[pivotRow]/pivotVal
//                M.columns[col].map { (row, Mij) -> negativeMagnitude(B[row] - Mij*k) }.sum()
//            }
//            .min()?:100
//    }


    // Returns -x*H(-x) where H is the Heaviside step function
    fun negativeMagnitude(x: Int) = if(x < 0) -x else 0



    // makes a pivot that brings the column 'col' closer to being
    // a "pivotable loop" which is a column with the properties:
    // - All non-zero elements are contained in X union Xb
    // - All elements are less than or equal to B
    // - There exists at least one element not in Xb that has value 1
    // Distance is measured as the L1 length of c-B + ABS(c-B) + |c intersect not X union Xb|
    //
    // Pivoting on an element with value +1 in a column other than 'col' has the effect of subtracting
    // r times the pivoted column to 'col', where r is the value in the same row in col (leaving the pivoted
    // row unchanged).
    //
    // Pivoting on an element with value -1 in a column other than 'col' has the effect of adding
    // r times the pivoted column to 'col' and flipping the sign of the pivoted row.
    //
    // OR...if we have a global measure of distance over all unpivoted Xb:
    // If we can pivot on a row in not X union Xb then do that
    // if we can pivot on a row in X not Xb, do that
    // otherwise we need to pivot in one element of Xb and pivot out another element in Xb
    // so choose one that reduces the overall distance
//    fun OldSimplex.pivotTowardsLoop(Xb: IntArrayVector, col: Int): Boolean {
//        return true
//    }


//    fun Simplex.pivotAtRandomValidRow(rows: Iterable<Int>): Boolean {
//        val validPivots = rows.flatMap { i ->
//            M.rows[i].keys
//                .map {j -> Simplex.PivotPoint(i,j)}
//                .filter { isValidPivotPoint(it) }
//        }
//        if(validPivots.isEmpty()) return false
//        val chosenPoint = validPivots.random()
//        println("Pivoting at $chosenPoint")
//        pivotAt(chosenPoint)
//        return true
//    }




}