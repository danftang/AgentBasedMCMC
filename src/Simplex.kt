import lib.sparseMatrix.*
import org.apache.commons.math3.distribution.EnumeratedDistribution
import org.apache.commons.math3.util.Pair
import java.lang.RuntimeException
import java.lang.StringBuilder
import java.util.AbstractMap
import kotlin.collections.HashMap
import kotlin.math.*
import kotlin.random.Random

// Represents MX = B
// with X is unknown and has no more non-zero elements than there are rows in M
class Simplex {
    var M: HashRowColIntMatrix                   // The matrix
    val B: HashIntVector                        // The observations and constraint values
    val X: HashIntVector                        // The solution MX = B
    val pivotPoints = ArrayList<Int>()          // the first pivotPoints.size rows are pivoted at column index stored in this
    val eventProbabilities: Map<Int,Double>
    var perturbationDistribution: List<Perturbation>? = null

    constructor(abmMatrix: SparseIntMatrix, observations: SparseIntVector, eventProbabilities: Map<Int,Double> = emptyMap()) {
        M = HashRowColIntMatrix(abmMatrix)
        B = HashIntVector(observations)
        X = HashIntVector()
        this.eventProbabilities = eventProbabilities
    }



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
    fun pivotTowards(Xb: HashIntVector): Boolean {
        val pivotedIn = pivotPoints.toSet()
        val colsToBePivotedIn = Xb.keys.filter { !pivotedIn.contains(it) }
        if(colsToBePivotedIn.isEmpty()) return false  // Xb already fully pivoted in
        // look for zero pivots
        colsToBePivotedIn.forEach { j ->
            M.columns[j]
                .find { (i, Mij) ->
                    Mij.absoluteValue == 1 && B[i] == 0 && Xb[pivotPoints[i]] == 0
                }
                ?.also { (i, _) ->
                    println("pivoting on zero")
                    pivotAt(PivotPoint(i,j))
                    return@pivotTowards true
                }
        }
        // look for loop pivots
        colsToBePivotedIn.forEach { j ->
            if(M.columns[j].all { (i, _) -> B[i] != 0 || Xb[pivotPoints[i]] != 0 }) {
                println("Found loop ${M.columns[j]} Xa members ${M.columns[j].filter { B[it.key] > 0 }} B = ${B.filter {M.columns[j].keys.contains(it.key)} }")
                M.columns[j]
                    .find { (i, Mij) -> Xb[pivotPoints[i]] == 0 && isValidPivotPoint(PivotPoint(i,j)) }
                    ?.also { (i, _) ->
                        println("pivoting on loop")
                        pivotAt(PivotPoint(i,j))
                        return@pivotTowards true
                    }
//                    ?:throw(RuntimeException("Found loop with no valid pivot point ${M.columns[j]}"))
            }
        }
        println("Cols to be pivoted in ${colsToBePivotedIn}")
        throw(RuntimeException("Can't find a zero pivot or a loop pivot. This shouldn't happen"))
    }

    // adds a false root and pivots to a single tree
    // Any pivot points in initialPivots are ensured to be
    // in the final tree
    fun pivotRootedSourceTree(initialPivots: Map<Int, Int> = emptyMap()) {
        val allSources = B.mapNotNull { if(it.value < 0) it.key else null }
        val rootRow = addSourceRoot(allSources)
        TreeConstructor(listOf(rootRow), initialPivots)
    }


    // adds a false root and pivots to a single tree
    // Any pivot points in initialPivots are ensured to be
    // in the final tree
    fun pivotSourcePolyTree(initialPivots: Map<Int, Int> = emptyMap()) {
        val allSources = B.mapNotNull { if(it.value < 0) it.key else null }
        TreeConstructor(allSources, initialPivots)
    }



    fun setPivotedSolutionValues(setNonPivotedValuesToZero: Boolean = true) {
        if(setNonPivotedValuesToZero) X.clear()
        for(pivotRow in pivotPoints.indices) {
            X[pivotPoints[pivotRow]] = B[pivotRow]
        }
    }


    fun pivotAt(pivot: PivotPoint) {
        applyPivot(pivot)
        if(pivot.row < pivotPoints.size) { // replacing an existing pivot
            pivotPoints[pivot.row] = pivot.col
        } else { // add newly pivoted row to the block for pivoted rows
            if (pivot.row != pivotPoints.size) swapRows(pivot.row, pivotPoints.size)
            pivotPoints.add(pivot.col)
        }
    }


    fun pivotAt(pivots: List<PivotPoint>) {
        pivots.forEach { pivot -> applyPivot(pivot) }

        val pivotsByRow = IntArray(M.nRows) { -1 }
        pivots.forEach { pivotsByRow[it.row] = it.col }
        for(i in pivotsByRow.indices) {
            if(pivotsByRow[i] != -1) {
                if(i < pivotPoints.size) { // replacing an existing pivot
                    pivotPoints[i] = pivotsByRow[i]
                } else { // add newly pivoted row to the block for pivoted rows
                    if (i != pivotPoints.size) swapRows(i, pivotPoints.size)
                    pivotPoints.add(pivotsByRow[i])
                }
            }
        }

//        var nextPivotRow = 0
//        while(nextPivotRow < M.nRows) {
//            while(nextPivotRow < M.nRows && pivotsByRow[nextPivotRow] == -1) nextPivotRow++
//            if(nextPivotRow < M.nRows) {
//                swapRows(pivotPoints.size, nextPivotRow)
//                pivotPoints.add(pivotsByRow[nextPivotRow])
//                nextPivotRow++
//            }
//        }
    }
//
//

    private fun applyPivot(pivot: PivotPoint) {
        if (M[pivot.row, pivot.col] == -1) {
            M.timesAssignRow(pivot.row, -1)
            B[pivot.row] *= -1
        }
        assert(M[pivot.row, pivot.col] == 1)
        val Bpiv = B[pivot.row]
        val Xpiv = if (pivot.row < pivotPoints.size) { // update pivoted X values
            val pivotOutCol = pivotPoints[pivot.row]
            val pivotOutVal = X[pivotOutCol]*M[pivot.row,pivotOutCol]
            X[pivot.col] = pivotOutVal
            X[pivotOutCol] = 0
            pivotOutVal
        } else {
            // while we're pivoting in, assume there are no semi-pivots
            X[pivot.col] = Bpiv
            Bpiv
        }

        assert(Bpiv == Xpiv)

        val col = M.columns[pivot.col].filter { it.key != pivot.row } // separate to prevent concurrent modification
        for (colEntry in col) {
            M.weightedRowPlusAssign(colEntry.key, pivot.row, -colEntry.value)
            if (Xpiv != 0 && colEntry.key < pivotPoints.size) { // might be distinct from B if non-pivoted values aren't zero
                X[pivotPoints[colEntry.key]] -= colEntry.value * Xpiv
            }
            B[colEntry.key] -= colEntry.value * Bpiv
        }
    }


    // perturbs the value of the unpivoted X elements, and updates the value
    // of the pivoted X elements accordingly
    fun columnPerturb() {
        val perturbationDistribution = findAllPositivePivotCols().map {
            Pair(it, 1.0)//exp(-changeInObservationNorm(it.value).toDouble()))
        }
        println("Found ${perturbationDistribution.size} possible perturbations")
        val chosenCol = EnumeratedDistribution(perturbationDistribution).sample()
        // update X_pivot for new X_unpivot value
        semiPivot(chosenCol, 1)
    }


    // perturbs the value of the unpivoted X element in the given row (col of M) by the given perturbation
    // and updates the value of the pivoted X elements accordingly
    fun semiPivot(j: Int, perturbation: Int) {
        X[j] += perturbation
        M.columns[j].forEach { X[pivotPoints[it.key]] -= it.value*perturbation }
    }


    fun pivotPerturb() {
        val possiblePivots = findAllValidPivotPoints().map {
//            Pair(it, 10.0.pow(-changeInObservationNorm(M[it.col]).toDouble()))
            Pair(it, if(changeInObservationNorm(M.columns[it.col]) == 0) 1.0 else 0.0)
        }
        println("Found ${possiblePivots.size} possible pivot point. Prob sum is ${possiblePivots.map {it.second}.sum()}")
        val choice = EnumeratedDistribution(possiblePivots).sample()
        pivotAt(choice)
        // setPivotedSolutionValues()
    }


    fun hybridPerturb() {
        val xPerturbations = positivePerturbations()
        val perturbationDistribution = (0 until M.nCols).asSequence()
            .flatMap { j ->
                xPerturbations[j]
                    ?.asSequence()
                    ?.mapNotNull { dXj -> if(dXj != 0) Pair(AbstractMap.SimpleEntry(j, dXj), 1.0) else null }
                    ?:emptySequence()
            }
            .toList()
        println("Found ${perturbationDistribution.size} possible perturbations")
        val chosenPerturbation = EnumeratedDistribution(perturbationDistribution).sample()
        perturbCol(chosenPerturbation.key, chosenPerturbation.value)
    }


    fun mcmcPerturb() {
        val perturbationDistribution = this.perturbationDistribution?:mcmcTransitionProbs()
        println("Found ${perturbationDistribution.size} possible perturbations")
        val chosenPerturbation = perturbationDistribution.sample()
        perturbCol(chosenPerturbation.column, chosenPerturbation.dx)
        // choose whether to accept perturbation
        val newPerturbationDistribution = mcmcTransitionProbs()
        val acceptance =  perturbationDistribution.sumOfProbabilities() / newPerturbationDistribution.sumOfProbabilities()
        if(acceptance < Random.nextDouble()) { // reject
            perturbCol(chosenPerturbation.column, -chosenPerturbation.dx)
            this.perturbationDistribution = perturbationDistribution
        } else {
            this.perturbationDistribution = newPerturbationDistribution
        }
    }


    // Minimise CX
    // Subject to:
    // MX = B
    // and
    // X >= 0 and X is integer
    // using Google OR-Tools
    fun MPsolve(C: List<Double> = DoubleArray(M.nCols) {0.0}.asList()) {
        X.setEqual(M.IPsolve(B, C, "=="))
    }


    class Perturbation(val column: Int, val dx: Int, val probability: Double, val cumulativeProb: Double)
    // returns the probability of perturbations of each unpivoted column,
    // given the current values of all other semi-pivoted columns
    fun mcmcTransitionProbs(): List<Perturbation> {
        var totalProb = 0.0
        val perturbations = (0 until M.nCols).asSequence()
            .flatMap { j ->
                if (isPivoted(j)) emptySequence() else {
                    val colProb = columnProbability(j)
                    M.columns[j]
                        .fold(IntRange(-X[j], Int.MAX_VALUE)) { range, (i, Mij) ->
                            val limit = X[pivotPoints[i]] / Mij
                            if (Mij > 0 && limit < range.last) IntRange(range.start, limit)
                            else if (Mij < 0 && limit > range.first) IntRange(limit, range.last)
                            else range
                        }
                        .asSequence()
                        .mapNotNull { dx ->
                            if(dx != 0) {
                                val prob = colProb.pow(0.5*dx) // square root to ensure acceptance equals ratio of sums
                                totalProb += prob
                                Perturbation(j, dx, prob, totalProb)
                            } else null
                        }
                }
            }.toList()
//        perturbations.forEach { perturbation -> // normalise probabilities
//            perturbation.probability /= totalProb
//            perturbation.cumulativeProb /= totalProb
//        }
        return perturbations
    }

    fun List<Perturbation>.sample(): Perturbation {
        val targetProb = Random.nextDouble(0.0, this.sumOfProbabilities())
        val index = binarySearch(0, size) { (it.cumulativeProb - targetProb).sign.toInt() }
        return this[index.absoluteValue]
    }

    fun List<Perturbation>.sumOfProbabilities(): Double = this.lastOrNull()?.cumulativeProb?:0.0

        // returns the probability of this column
    // i.e. the product of its event probabilities raised to the power of their occupation number
    fun columnProbability(j: Int): Double {
        return M.columns[j].fold(1.0) { prob, entry ->
            prob * (eventProbabilities[entry.key]?.pow(entry.value)?:0.0)
        }
    }


    // Sets column j to the given value,
    // doing a full pivot where possible or semi-pivot otherwise
    fun perturbCol(j: Int, perturbation: Int) {
        assert(!isPivoted(j))
        M.columns[j].find { (i, Mij) ->  // try to find a possible pivot point
            Mij.absoluteValue == 1 && i < pivotPoints.size && X[pivotPoints[i]] == perturbation
        }?.run {
            pivotAt(PivotPoint(key, j))
        }?:run {
            semiPivot(j, perturbation)
        }
    }


    // all state-changing pivot points that will not produce negative act values
    fun findAllValidPivotPoints(): List<PivotPoint> {
        return B.flatMap { Bentry ->
            if(Bentry.key < pivotPoints.size && Bentry.value > 0) {
                M.rows[Bentry.key].mapNotNull { MTentry ->
                    val pivot = PivotPoint(Bentry.key, MTentry.key)
                    if (isValidPivotPoint(pivot)) pivot else null
                }
            } else emptyList()
        }
    }


//    // finds all pairs of columns that when added together
//    //
//    fun findAllPerturbationColumnPairs() {
//
//    }

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


//    fun findInitialSolution() {
//        val initialSolution = M.IPsolve(B, constraintType = "=")
//        val Xinit = SparseColIntMatrix.SparseIntColumn(initialSolution)
//        assert(M*Xinit == B)
//        Xinit.keys.forEach { M[it].keys.forEach { assert(it < pivotPoints.size) } }
//        for(j in initialSolution.indices) {
//            if(initialSolution[j] != 0) {
//                if(M[j].sparseSize > 1) {
//                    println("Pivoting...")
//                    val pivot = PivotPoint(
//                        M[j].entries.find { it.value == 1 && it.key < pivotPoints.size }?.key?:throw(RuntimeException("Can't pivot solution. Odd!"))
//                        ,j
//                    )
//                    applyPivot(pivot)
//                    pivotPoints[pivot.row] = pivot.col
//                    M = MT.transpose()
//                    println("Observation error ${observationNorm()}")
//                }
//            }
//        }
//        setPivotedSolutionValues()
//        // TODO: just a sanity check:
//        assert(M*X == B)
//    }


    // pivots in all non-zero columns in the supplied solution
    // then greedily pivots in any remaining rows
    fun pivotToSolution(solution: SparseIntVector) {
        assert(M*solution == B)
//        solution.keys.forEach { act -> M[act].keys.forEach { row ->
//            println("$row OK")
//            assert(row < pivotPoints.size)
//        } }
        for(entry in solution) {
            if(M.columns[entry.key].sparseSize > 1) {
//                println("Pivoting...")
                val pivot = PivotPoint(
                    M.columns[entry.key].find { it.value.absoluteValue == 1 && (it.key >= pivotPoints.size || solution[pivotPoints[it.key]] == 0) }?.key?:throw(RuntimeException("Can't pivot to solution. Odd!"))
                    ,entry.key
                )
                pivotAt(pivot)
//                println("Observation error ${observationNorm()}")
            }
        }
        greedyPivotAll()
        setPivotedSolutionValues(true)
    }


    // pivots in any rows that aren't already pivoted in
    // using a greedy algorithm:
    // for each row in turn, choose to pivot on the column with smallest
    // sparse size among all that have a +-1 in this row
    fun greedyPivotAll() {
        for(pivotRow in pivotPoints.size until M.nRows) {
            val pivotCol = M.rows[pivotRow]
                .asSequence()
                .filter { it.value.absoluteValue == 1 }
                .minBy { M.columns[it.key].sparseSize }
                ?.run { this.key }
                ?:throw(RuntimeException("No unit elements to pivot on"))
            pivotAt(PivotPoint(pivotRow,pivotCol))
        }
    }


    fun pivotOutNegatives() {
        var tolerance = 0
        while(negativeActCount() > 0) {
            println("negative score is ${negativeActCount()}")
            allPivotPoints()
                .find { changeInNegativeActCountOnPivot(it) < tolerance }
                ?.run {
                    println("Pivoting out negative at $this")
                    val expectedChange = changeInNegativeActCountOnPivot(this)
                    println("expected count change is $expectedChange")
                    val startScore = negativeActCount()
                    assert(M*X == B)
                    pivotAt(this)
                    println("actual change is ${negativeActCount() - startScore}")
                    assert(M*X == B)
                    assert(negativeActCount() - startScore == expectedChange )
                    tolerance = 0
                }
                ?:run {
                    println("increasing tolerance")
                    tolerance++
                    if(tolerance == 3) throw(RuntimeException("No possible pivots to pivot out negatives"))
                }
        }
        // setPivotedSolutionValues()
    }

    fun negativeActCount() = X.values.sumBy { max(0,-it) }
        //B.sumBy { if(it.key < pivotPoints.size && it.value < 0) 1 else 0  }


    // calculates the change in negativeActCount upon pivot at the given point
    fun changeInNegativeActCountOnPivot(pivot: PivotPoint): Int {
        val Bpiv = B[pivot.row]
    //    assert(Bpiv >= 0)
        val Mpiv = M[pivot.row, pivot.col]
        assert(Mpiv.absoluteValue == 1)
        return M.columns[pivot.col].sumBy { entry ->
            if(entry.key != pivot.row && entry.key < pivotPoints.size) {
                val bi = B[entry.key]
                max(0, Bpiv*entry.value/Mpiv - bi) - max(0, -bi)
            } else 0
        } + if(Mpiv < 0) Bpiv else 0
    }






    // returns all points that have absolute value one and B value not equal to zero
    fun allPivotPoints(): Sequence<PivotPoint> {
        return (0 until pivotPoints.size).asSequence()
            .flatMap { row ->
                if(B[row] != 0) {
                    M.rows[row].asSequence().mapNotNull { (col, Mij) ->
                        if(Mij.absoluteValue == 1 && col != pivotPoints[row]) PivotPoint(row, col) else null
                    }
                } else emptySequence()
            }
    }

    // finds pivot cols that are non-zero in given row
    fun findPositivePivotCols(row: Int) =
        M.columns
            .withIndex()
            .filter {
                it.value.sparseSize > 1 && it.value.keys.contains(row) && isPositiveColumn(it.value) //&& changeInObservationNorm(it.value) <= 0
            }


    // returns the range of possible perturbations of each unpivoted column,
    // given the current values of all other semi-pivoted columns
    // or null if the column is already pivoted
    fun positivePerturbations(): List<IntRange?> {
        return Array(M.nCols) { j ->
            if(isPivoted(j)) null else {
                M.columns[j].asSequence()
                    .fold(IntRange(-X[j], Int.MAX_VALUE)) { range, (i, Mij) ->
                        val limit = X[pivotPoints[i]]/Mij
                        if(Mij > 0 && limit < range.last) IntRange(range.start, limit)
                        else if(Mij < 0 && limit > range.first) IntRange(limit, range.last)
                        else range
                    }
            }
        }.asList()
    }




    // true if this column is currently pivoted
    fun isPivoted(j: Int): Boolean {
        val col = M.columns[j]
        if(col.sparseSize != 1) return false
        val i = col.first().key
        return i < pivotPoints.size && pivotPoints[i] == j
    }


    fun findAllPositivePivotCols(): List<Int> =
        (0 until M.nCols).filter {
            val col = M.columns[it]
            col.sparseSize > 1 && isPositiveColumn(col)
        }
//        M.columns.asSequence()
//            .withIndex()
//            .filter {
//                it.value.sparseSize > 1 && isPositiveColumn(it.value) //&& changeInObservationNorm(it.value) <= 0
//            }
//            .toList()


    // how desirable is the perturbation of X[col] to X[col] + perturbation
    // where col >= pivotpoints.size
    //
    fun perturbationScore(col: Int, perturbation: Int): Double {
        return 1.0
//        var score = 0
//        for(entry in M.columns[col]) {
//            if(entry.key < pivotPoints.size) {
//                if(pivotPoints[entry.key])
//                if(X[pivotPoints[entry.key]] < entry.value) score++
//
//            } else -it.key
//
//        }
//        val Xperturbation = M.columns[col].mapKeys {
//            if(it.key < pivotPoints.size) pivotPoints[it.key] else -it.key
//        }
//        return score(X - Xperturbation)
    }


    fun score(x: SparseIntVector): Double {
        val rootOffset = M.rows[M.nRows-1].dotProd(x).absoluteValue
        return -(x.count {it.key != 1-M.nRows && (it.key < 0 || it.value < 0) } + rootOffset).toDouble()
    }

    // true iff subtracting this column from B increases the 1-norm of the observation
    fun changeInObservationNorm(col: SparseIntVector): Int {
        var dErr = 0
        for(entry in col) {
            if(entry.key >= pivotPoints.size) {
                val bi = B[entry.key]
                dErr += (bi - entry.value).absoluteValue - bi.absoluteValue
            }
        }
        return dErr
    }

    // returns the 1-norm of the observation
//    fun observationNorm(): Int {
//        var norm = 0
//        for(entry in B) {
//            if(entry.key >= pivotPoints.size) {
//                norm += entry.value.absoluteValue
//            }
//        }
//        return norm
//    }



    // true iff setting this column to one will not create any negative values in B
    fun isPositiveColumn(col: SparseIntVector): Boolean {
        return col.all { entry ->
//                B[entry.key] >= entry.value
                entry.key < pivotPoints.size && X[pivotPoints[entry.key]] >= entry.value
            }
    }

    // true iff this point equals +-1, is not already pivoted-in and
    // pivoting at this point will not create any negative values in B
    fun isValidPivotPoint(pivot: PivotPoint): Boolean {
        val pivCol = M.columns[pivot.col]
        if(pivCol[pivot.row].absoluteValue != 1 || pivCol.sparseSize <= 1) return false
        val k = B[pivot.row]/pivCol[pivot.row]
        return k >= 0 && pivCol.all { B[it.key] >= k*it.value || it.key >= pivotPoints.size }
    }


    fun swapRows(i1: Int, i2: Int) {
        if((i1 < pivotPoints.size) xor (i2 < pivotPoints.size)) throw(RuntimeException("Can't swap pivoted row with non-pivoted row"))
        if(i1 == i2) return
        M.swapRows(i1,i2)
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
    fun toCompositeMatrix(): HashRowColIntMatrix {
        val compositeMatrix = M.copy()
        compositeMatrix.resize(M.rows.size, M.columns.size + 1)
        compositeMatrix.setColumn(M.columns.size, B)
        return compositeMatrix
    }

    override fun toString(): String {
        val string = StringBuilder()
        val activeCols = CharArray(M.nCols*3) { '.' }
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
        val activeCols = CharArray(M.nCols) { '.' }
        pivotPoints.forEach { j ->
            activeCols[j] = '#'
        }
        string.append(activeCols)
        string.append('\n')
        string.append(toCompositeMatrix().toSparsityString())
        return string.toString()

    }


    // Adds a 'fake root node that has value zero and edges from it to each
    // source node, so that the whole matrix becomes a single tree with only forward edges
    // takes a list of row indices which identify the source rows
    // and adds the fake root to the bottom of this matrix (increasing nRows by 1)
    // and a new action for each source (increasing nCols by the number of sources)
    // returns the row index of the added root
    fun addSourceRoot(sources: List<Int>): Int {
        val i = M.nRows
        val j = M.nCols
        M.resize(i + 1, j + sources.size)
        for(col in sources.indices) {
            M[i, j + col] = -1
            M[sources[col], j + col] = 1
        }
        return M.nRows-1
    }


    data class PivotPoint(val row: Int, val col: Int)

//    enum class RowType {
//        PIVOTED,
//        PIVOTEDROOT,
//        UNPIVOTED,
//        UNPIVOTEDROOT
//
//        fun isPivoted() = (this == PIVOTED || this == PIVOTEDROOT)
//    }

    inner class TreeConstructor {
        val ROOT = -2
        val distanceUpperBounds =
            HashMap<Int, MutableList<PivotPoint>>() // upper bounds of distance to root. Maps distance to pivot points
        val colDistances: IntArray              // lower bound on distance when pivoting on this column
                                                // when the number of undecided rows for a column becomes 1, this bound is tight
        var distanceLowerBound = 0
        var highestBound = 0

        constructor(rootRows: List<Int>, initialPivots: Map<Int,Int> = emptyMap()) {
            distanceUpperBounds[0] = rootRows.map { PivotPoint(it, ROOT) }.toMutableList()

            val Mcols = HashColIntMatrix(M)
            colDistances = IntArray(M.nCols) { 1 }
            while(distanceLowerBound <= highestBound) {
                val leaves = distanceUpperBounds[distanceLowerBound]
                    ?.filter { it.col == ROOT || Mcols[it.row, it.col] != 0 } // remove if already reduced
                    ?.distinctBy { it.row }
                    ?:listOf()
                distanceUpperBounds[distanceLowerBound] = leaves.toMutableList()
                leaves
                    .forEach { pivotPoint -> // recalculate colDistances and add new children to upper bounds
                        M.rows[pivotPoint.row]
                            .forEach { (childCol, childVal) ->
                                assert(Mcols[pivotPoint.row, childCol] != 0)
                                val priorSize = Mcols.columns[childCol].sparseSize
//                                println("initial val = ${M[pivotPoint.row, childCol]}")
                                Mcols[pivotPoint.row, childCol] = 0 // use Mcols to keep track of followed values
//                                println("prior size = $priorSize new size = ${M.columns[childCol].sparseSize}")
                                assert(Mcols.columns[childCol].sparseSize == priorSize - 1)
                                colDistances[childCol] += distanceLowerBound
                                if(Mcols.columns[childCol].sparseSize == 1 && childVal < 0) { // follow children forward in time when reduced
                                    distanceUpperBounds
                                        .getOrPut(colDistances[childCol]) {ArrayList()}
                                        .add(PivotPoint(Mcols.columns[childCol].keys.first(), childCol))
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
            greedyPivotAll() // finally, pivot on tree roots
        }

    }
}