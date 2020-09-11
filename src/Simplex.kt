import lib.SparseColIntMatrix
import org.apache.commons.math3.distribution.EnumeratedDistribution
import org.apache.commons.math3.util.Pair
import java.lang.RuntimeException
import java.lang.StringBuilder
import kotlin.collections.HashMap
import kotlin.math.*

// Represents MX = B
// with X is unknown and has no more non-zero elements than there are rows in M
class Simplex {
    var MT: SparseColIntMatrix                  // The matrix, transposed
    var M: SparseColIntMatrix                   // The matrix, (not transposed)
    val B: SparseColIntMatrix.SparseIntColumn   // The current solution values
    val X: SparseColIntMatrix.SparseIntColumn   // The solution
    val pivotPoints = ArrayList<Int>()          // first pivotPoints.size rows are pivoted at column index stored in this

    val nRows: Int
        get() = MT.nCols

    val nCols: Int
        get() = M.nCols

    constructor(abmMatrix: SparseColIntMatrix, observations: SparseColIntMatrix.SparseIntColumn) {
        M = abmMatrix.copy()
        MT = abmMatrix.transpose()
        B = observations.copy()
        X = SparseColIntMatrix.SparseIntColumn()
    }

    fun pivotPolyTree(initialPivots: Map<Int, Int> = emptyMap()) {
        TreeConstructor(
            B.entries
                .filter { it.value < 0 }    // only sources for roots
                .map { it.key }
            , initialPivots
        )
    }


    fun setPivotedSolutionValues(setNonPivotedValuesToZero: Boolean = true) {
        if(setNonPivotedValuesToZero) X.data.clear()
        for(pivotRow in pivotPoints.indices) {
            X[pivotPoints[pivotRow]] = B[pivotRow]
        }
    }


    fun pivotAt(pivot: PivotPoint) { pivotAt(listOf(pivot)) }


    fun pivotAt(pivots: List<PivotPoint>) {
        pivots.forEach { pivot ->
            applyPivot(pivot)
        }
        val pivotsByRow = IntArray(nRows) { -1 }
        pivots.forEach { pivotsByRow[it.row] = it.col }
        var nextPivotRow = 0
        while(nextPivotRow < nRows) {
            while(nextPivotRow < nRows && pivotsByRow[nextPivotRow] == -1) nextPivotRow++
            if(nextPivotRow < nRows) {
                swapRows(pivotPoints.size, nextPivotRow)
                pivotPoints.add(pivotsByRow[nextPivotRow])
                nextPivotRow++
            }
        }
    }


    private fun applyPivot(pivot: PivotPoint) {
        val pivotEq = MT[pivot.row]
//        if(pivotEq[pivot.col] == -1) {
//            pivotEq *= -1
//            B[pivot.row] *= -1
//        }
        assert(pivotEq[pivot.col] == 1)
        val Bpiv = B[pivot.row]
        for(i in 0 until nRows) {
            if(i != pivot.row) {
                val selectedColVal = MT[i][pivot.col]
                if (selectedColVal != 0) {
                    MT[i].weightedPlusAssign(pivotEq, -selectedColVal)
                    B[i] -= selectedColVal * Bpiv
                }
            }
        }
    }



//    fun perturb() {
//        val perturbationDistribution = findAllPositivePivotCols().map {
//            Pair(it, exp(-changeInObservationNorm(it.value).toDouble()))
//        }
////        println("possible perturbations $perturbationDistribution")
//        println("observation norm perturbations ${perturbationDistribution.map { -ln(it.second).toInt() }}")
//        val choice = EnumeratedDistribution(perturbationDistribution).sample()
////        println("perturbing with col ${choice.value}")
////        println("initial B = $B")
//        X[choice.index] += 1
//        B -= choice.value
////        println("new B = $B")
//        setPivotedSolutionValues(false)
//        println("new observation norm is ${observationNorm()}")
//    }

    fun pivotPerturb() {
        val possiblePivots = findAllValidPivotPoints().map {
//            Pair(it, 10.0.pow(-changeInObservationNorm(M[it.col]).toDouble()))
            Pair(it, if(changeInObservationNorm(M[it.col]) == 0) 1.0 else 0.0)
        }
        println("Found ${possiblePivots.size} possible pivot point. Prob sum is ${possiblePivots.map {it.second}.sum()}")
        val choice = EnumeratedDistribution(possiblePivots).sample()
        applyPivot(choice)
        M = MT.transpose()
        pivotPoints[choice.row] = choice.col
        setPivotedSolutionValues()
    }

    // all pivot points that will not produce negative act values
    fun findAllValidPivotPoints(): List<PivotPoint> {
        return B.entries.flatMap { Bentry ->
            if(Bentry.key < pivotPoints.size && Bentry.value > 0) {
                MT[Bentry.key].entries.mapNotNull { MTentry ->
                    val pivot = PivotPoint(Bentry.key, MTentry.key)
                    if (isValidPivotPoint(pivot)) pivot else null
                }
            } else emptyList()
        }
    }


//    fun removeNegatives() {
//            while(B.entries.find { it.value < 0 }?.apply {
//                    val perturbations = findPositivePivotCols(key)
//                        .sortedBy { changeInObservationNorm(it.value) }
//                    val chosenCol = perturbations[0]
//                    X[chosenCol.index] += 1
//                    B -= chosenCol.value
//                    setPivotedSolutionValues(false)
//                    println("added col ${chosenCol.value}")
//                    println("observation norm is ${observationNorm()}")
//                    println("new solution is $X")
//            } != null) {}
//    }


    fun findInitialSolution() {
        M = MT.transpose()
        val initialSolution = M.IPsolve(B, constraintType = "=")
        val Xinit = SparseColIntMatrix.SparseIntColumn(initialSolution)
        assert(M*Xinit == B)
        Xinit.keys.forEach { M[it].keys.forEach { assert(it < pivotPoints.size) } }
        for(j in initialSolution.indices) {
            if(initialSolution[j] != 0) {
                if(M[j].sparseSize > 1) {
                    println("Pivoting...")
                    val pivot = PivotPoint(
                        M[j].entries.find { it.value == 1 && it.key < pivotPoints.size }?.key?:throw(RuntimeException("Can't pivot solution. Odd!"))
                        ,j
                    )
                    applyPivot(pivot)
                    pivotPoints[pivot.row] = pivot.col
                    M = MT.transpose()
                    println("Observation error ${observationNorm()}")
                }
            }
        }
        setPivotedSolutionValues()
        // TODO: just a sanity check:
        assert(M*X == B)
    }


    fun pivotToSolution(solution: SparseColIntMatrix.SparseIntColumn) {
        assert(M*solution == B)
//        solution.keys.forEach { act -> M[act].keys.forEach { row ->
//            println("$row OK")
//            assert(row < pivotPoints.size)
//        } }
        for(entry in solution.entries) {
            if(M[entry.key].sparseSize > 1) {
//                println("Pivoting...")
                val pivot = PivotPoint(
                    M[entry.key].entries.find { it.value == 1 && it.key < pivotPoints.size }?.key?:throw(RuntimeException("Can't pivot to solution. Odd!"))
                    ,entry.key
                )
                applyPivot(pivot)
                pivotPoints[pivot.row] = pivot.col
                M = MT.transpose()
//                println("Observation error ${observationNorm()}")
            }
        }
        setPivotedSolutionValues()
    }


    fun pivotOutNegatives() {
        var tolerance = 0
        while(negativeActCount() > 0) {
            println("negative score is ${negativeActCount()}")
            allPivotPoints()
                .find { changeInNegativeActCountOnPivot(it) < tolerance }
                ?.run {
                    println("Pivoting out negative at $this")
                    println("expected count change is ${changeInNegativeActCountOnPivot(this)}")
                    println("negative score is ${negativeActCount()}")
                    applyPivot(this)
                    M = MT.transpose()
                    pivotPoints[this.row] = this.col
                    println("new negative score is ${negativeActCount()}")
                    tolerance = 0
                }
                ?:run {
                    println("increasing tolerance")
                    tolerance++
                    if(tolerance == 3) throw(RuntimeException("No possible pivots to pivot out negatives"))
                }
        }
        setPivotedSolutionValues()
    }

    fun negativeActCount() = B.entries.sumBy { if(it.key < pivotPoints.size && it.value < 0) 1 else 0  }

    // calculates the change in negativeActCount upon pivot at the given point
    fun changeInNegativeActCountOnPivot(pivot: PivotPoint): Int {
        val Bpiv = B[pivot.row]
        assert(Bpiv >= 0)
        return M[pivot.col].entries.sumBy { entry ->
            if(entry.key != pivot.row && entry.key < pivotPoints.size) {
                val bi = B[entry.key]
                max(0, Bpiv*entry.value - bi) - max(0, -bi)
            } else 0
        }
    }


    // returns all points that have value one and are on a row with positive B
    fun allPivotPoints(): Sequence<PivotPoint> {
        return (0 until pivotPoints.size).asSequence()
            .flatMap { row ->
                if(B[row] > 0) {
                    MT[row].entries.asSequence().mapNotNull { MTentry ->
                        if(MTentry.value == 1 && MTentry.key != pivotPoints[row]) PivotPoint(row, MTentry.key) else null
                    }
                } else emptySequence()
            }
    }

    // finds pivot cols that are non-zero in given row
    fun findPositivePivotCols(row: Int) =
        M
            .withIndex()
            .filter {
                it.value.sparseSize > 1 && it.value.keys.contains(row) && isPositiveColumn(it.value) //&& changeInObservationNorm(it.value) <= 0
            }


    fun findAllPositivePivotCols(): List<IndexedValue<SparseColIntMatrix.SparseIntColumn>> =
        M
            .withIndex()
            .filter {
                it.value.sparseSize > 1 && isPositiveColumn(it.value) //&& changeInObservationNorm(it.value) <= 0
            }


    // true iff subtracting this column from B increases the 1-norm of the observation
    fun changeInObservationNorm(col: SparseColIntMatrix.SparseIntColumn): Int {
        var dErr = 0
        for(entry in col.entries) {
            if(entry.key >= pivotPoints.size) {
                val bi = B[entry.key]
                dErr += (bi - entry.value).absoluteValue - bi.absoluteValue
            }
        }
        return dErr
    }

    // returns the 1-norm of the observation
    fun observationNorm(): Int {
        var norm = 0
        for(entry in B.entries) {
            if(entry.key >= pivotPoints.size) {
                norm += entry.value.absoluteValue
            }
        }
        return norm
    }



    // true iff setting this column to one will not create any negative values in B
    fun isPositiveColumn(col: SparseColIntMatrix.SparseIntColumn): Boolean {
        return col.entries
            .fold(true) { isPositive, entry ->
                isPositive && B[entry.key] >= entry.value
            }
    }


    // true iff pivoting at this point will not create any negative values in B
    fun isValidPivotPoint(pivot: PivotPoint): Boolean {
        val pivCol = M[pivot.col]
        if(pivCol[pivot.row] != 1 || pivCol.sparseSize <= 1) return false
        val Bpiv = B[pivot.row]
        return pivCol.entries
            .fold(true) { isPositive, entry ->
                isPositive && (B[entry.key] >= Bpiv*entry.value || entry.key >= pivotPoints.size)
            }
    }


    fun swapRows(i1: Int, i2: Int) {
        if((i1 < pivotPoints.size) xor (i2 < pivotPoints.size)) throw(RuntimeException("Can't swap pivoted row with non-pivoted row"))
        if(i1 == i2) return
        MT.swapCols(i1,i2)
        val tmp = B[i1]
        B[i1] = B[i2]
        B[i2] = tmp
        if(i1 < pivotPoints.size) {
            val ppTmp = pivotPoints[i1]
            pivotPoints[i1] = pivotPoints[i2]
            pivotPoints[i2] = ppTmp
        }
    }

    // returns M with B added as a rightmost column
    fun toCompositeMatrix(): SparseColIntMatrix {
        val compositeMatrix = MT.transpose()
        compositeMatrix.addColRight(B)
        return compositeMatrix
    }

    override fun toString(): String {
        val string = StringBuilder()
        val activeCols = CharArray(nCols*3) { '.' }
        pivotPoints.forEach { j ->
            activeCols[j*3 + 1] = '#'
        }
        string.append(activeCols)
        string.append('\n')
        string.append(toCompositeMatrix().toString())
        return string.toString()
    }

    fun toSparsityString(): String {
        val string = StringBuilder()
        val activeCols = CharArray(nCols) { '.' }
        pivotPoints.forEach { j ->
            activeCols[j] = '#'
        }
        string.append(activeCols)
        string.append('\n')
        string.append(toCompositeMatrix().toSparsityString())
        return string.toString()

    }

    data class PivotPoint(val row: Int, val col: Int)

    inner class TreeConstructor {
        val ROOT = -2
        val distanceUpperBounds =
            HashMap<Int, MutableList<PivotPoint>>() // upper bounds of distance to root. Maps distance to pivot points
        val colDistances: IntArray              // lower bound on distance when pivoting on this column
                                                // when the number of undecided rows for a column becomes 1, this bound is tight
        var distanceLowerBound = 0
        var highestBound = 0

        constructor(rootRows: List<Int>, initialPivots: Map<Int,Int> = emptyMap()) {
            // add roots as zero distance upper bounds
            addSourceRoot(rootRows)
            distanceUpperBounds[0] = mutableListOf(PivotPoint(nRows-1, ROOT)) // rootRows.map { PivotPoint(it,ROOT) }.toMutableList()

            colDistances = IntArray(nCols) { 1 }
            while(distanceLowerBound <= highestBound) {
                val leaves = distanceUpperBounds[distanceLowerBound]
                    ?.filter { it.col == ROOT || M[it.row, it.col] != 0 } // remove if already reduced
                    ?.distinctBy { it.row }
                    ?:listOf()
                distanceUpperBounds[distanceLowerBound] = leaves.toMutableList()
                leaves
                    .forEach { pivotPoint -> // recalculate colDistances and add new children to upper bounds
                        MT[pivotPoint.row].entries
                            .forEach { (childCol, childVal) ->
                                assert(MT[childCol, pivotPoint.row] != 0)
                                assert(M[pivotPoint.row, childCol] != 0)
                                val priorSize = M[childCol].sparseSize
                                M[pivotPoint.row, childCol] = 0 // use M to keep track of followed values
                                assert(M[childCol].sparseSize == priorSize - 1)
                                colDistances[childCol] += distanceLowerBound
                                if(M[childCol].sparseSize == 1 && childVal < 0) { // follow children forward in time when reduced
                                    distanceUpperBounds
                                        .getOrPut(colDistances[childCol]) {ArrayList()}
                                        .add(PivotPoint(M[childCol].keys.first(), childCol))
                                    if(colDistances[childCol] > highestBound) highestBound = colDistances[childCol]
                                }
                            }
                    }
                distanceLowerBound++
            }

//            assert(highestBound <= distanceLowerBound)

//            println("Pivot map is $distanceUpperBounds")
            // now we can read off the tree from the distanceUpperBounds
            // pivoting in reverse order is most efficient

            val pivotsInOrder = (distanceLowerBound-1 downTo 1)
                .flatMap {
                    distanceUpperBounds[it]
                        ?.map { originalPivot ->
                            initialPivots[originalPivot.row]
                                ?.run { PivotPoint(originalPivot.row, this) }
                                ?:originalPivot
                        }
                        ?:emptyList()
                }

            pivotAt(pivotsInOrder)
            M = MT.transpose()
            setPivotedSolutionValues()
        }

        // Adds a 'fake root node that has value zero and edges from it to each
        // source node, so that the whole matrix becomes a single tree with only forward edges
        // takes a list of row indices which identify the source rows
        // and adds the fake root to the bottom of this matrix (increasing nRows by 1)
        // and a new action for each source (increasing nCols by the number of sources)
        fun addSourceRoot(sources: List<Int>) {
            for(sourceRow in sources) {
                M.addColRight(SparseColIntMatrix.SparseIntColumn(
                    sourceRow to 1,
                    M.nRows to -1
                ))
            }
            M.nRows += 1
            MT = M.transpose()
        }

    }
}