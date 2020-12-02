import lib.*
import lib.sparseIntMatrix.HashIntVector
import lib.sparseIntMatrix.HashRowColIntMatrix
import lib.sparseIntMatrix.SparseIntMatrix
import lib.sparseIntMatrix.SparseIntVector
import org.apache.commons.math3.fraction.Fraction
import org.apache.commons.math3.util.ArithmeticUtils.gcd
import org.apache.commons.math3.util.Pair
import java.lang.RuntimeException
import java.lang.StringBuilder
import kotlin.math.*
import kotlin.random.Random

// Represents MX = B
// with X is unknown and has no more non-zero elements than there are rows in M
//class OldSimplex {
//    var M: HashRowColIntMatrix                   // The matrix
////    val B: HashIntVector                        // The observations and constraint values
//    val X: Map<Int,Fraction>                        // The solution MX = B
//    val basicColsByRow: IntArray          // gives the column that is currently pivoted in on this row, or -1 if not yet pivoted
//    val eventProbabilities: Map<Int,Double>
//    var perturbationDistribution: List<Perturbation>? = null
//
//    constructor(abmMatrix: SparseIntMatrix, eventProbabilities: Map<Int,Double> = emptyMap(), initialSolution: SparseIntVector) {
//        M = HashRowColIntMatrix(abmMatrix)
//        X = HashMap<Int,Fraction>()
//        initialSolution.forEach {
//            X[it.key] = Fraction(it.value,1)
//        }
//        this.eventProbabilities = eventProbabilities
//        // Now set up initial pivot
//        basicColsByRow = IntArray(abmMatrix.nRows) { -1 }
//        pivotToSolution(initialSolution)
//    }
//
//
//    constructor(abmMatrix: SparseIntMatrix, observations: SparseIntVector, eventProbabilities: Map<Int,Double> = emptyMap()) : this(
//        abmMatrix,
//        eventProbabilities,
//        abmMatrix.IPsolve(observations, DoubleArray(abmMatrix.nCols) {0.0}.asList(), "==").toSparseIntVector()
//    )
//
//    // returns the range of perturbations that the "j"'th element of X
//    // can take while maintaining a positive solution
//    // given that all other non-basic elements remain fixed
//    // If X_j is a basic element, returns the null range 0..0
//    // otherwise, for a perturbation DX_j we have
//    // Mip[i]*dX_p[i] = -Mij*dX_j
//    // so
//    // dX_p[i]/dX_j = -Mij/Mip[i]
//    // so DX_j is at a limit for the i'th row when
//    // dX_p[i] = -X_p[i]
//    // and
//    // dX_j(i) = X_p[i]*Mip[i]/Mij
//    // since X_p[i] and Mip[i] are +ve, DX_j(i) has the same sign as Mij
//    // so the range of dXj is
//    // max(-X_j,arg max_(i s.t. Mij<0)(Mip[i]*X_p[i]/Mij)) ... arg min_(i s.t. Mij>0)(Mip[i]*X_p[i]/Mij)
//    fun perturbationRange(j: Int): Pair<Fraction,Fraction> {
//        if(isBasicColumn(j)) return Pair(Fraction.ZERO, Fraction.ZERO)
//        var dXjmin = Fraction(Int.MIN_VALUE, 1)
//        var dXjmax = Fraction(Int.MAX_VALUE,1)
//        M.columns[j].forEach { (i, Mij) ->
//            val pi = basicColsByRow[i]
//            val dXji = (X[pi]*M[i,pi])/Mij
//            if(Mij > 0) {
//                if (dXji < dXjmax) dXjmax = dXji
//            } else {
//                if(dXji > dXjmin) dXjmin = dXji
//            }
//        }
//        return Pair(max(-X[j], dXjmin), dXjmax)
//    }
//
//
//    fun maxPerturbation(j: Int): Fraction {
//        var dXjmax = Fraction(Int.MAX_VALUE,1)
//        M.columns[j].forEach { (i, Mij) ->
//            val pi = basicColsByRow[i]
//            val dXji = (X[pi]*M[i,pi])/Mij
//            if(Mij > 0 && dXji < dXjmax) dXjmax = dXji
//        }
//        return dXjmax
//    }
//
//
//    // Returns the rows in column j that are
//    // pivot points that maintain a positive solution
//    fun pivotableRows(j: Int): List<Int> {
//        if(isBasicColumn(j)) return emptyList()
//
//        var dXjmax = Fraction(Int.MAX_VALUE,1)
//        val limits = M.columns[j].map { (i, Mij) ->
//            val pi = basicColsByRow[i]
//            val dXji = (X[pi]*M[i,pi])/Mij
//            if(Mij > 0 && dXji < dXjmax) dXjmax = dXji
//            Pair(i,dXji)
//        }
//
//        return limits
//            .filter { it.second == dXjmax }
//            .map { it.first }
//    }
//
//    // Pivots-in column j. Choosing a column to pivot-out
//    // at random from among the columns that can be pivoted
//    // out while maintaining a positive solution. If no
//    // such column exists, then returns false
//    fun pivotUp(j: Int): Boolean {
//        return false
//    }
//
//    // adds a false root and pivots to a single tree
//    // Any pivot points in initialPivots are ensured to be
//    // in the final tree
////    fun pivotRootedSourceTree(initialPivots: Map<Int, Int> = emptyMap()) {
////        val allSources = B.mapNotNull { if(it.value < 0) it.key else null }
////        val rootRow = addSourceRoot(allSources)
////        TreeConstructor(listOf(rootRow), initialPivots)
////    }
//
//
//    // Pivots to a poly-tree where each source row (B value < 1)
//    // is a root
//    // Any pivot points in initialPivots are ensured to be
//    // in the final tree
////    fun pivotSourcePolyTree(initialPivots: Map<Int, Int> = emptyMap()) {
////        val allSources = B.mapNotNull { if(it.value < 0) it.key else null }
////        TreeConstructor(allSources, initialPivots)
////    }
//
//
//    // Sets the basic variables in X, and set all semi-pivots to zero
////    fun setXToBasicValue() {
////        X.clear()
////        for(pivotRow in basicColsByRow.indices) {
////            val pivotCol = basicColsByRow[pivotRow]
////            X[pivotCol] = Fraction(B[pivotRow],M[pivotRow,pivotCol])
////        }
////    }
//
//    // performs a pivot and updates X value if not an initial pivot-in
//    fun pivotAt(pivot: PivotPoint) {
//        val pivotOutCol = basicColsByRow[pivot.row]
//        applyPivot(pivot)
//        applyNullVector(pivotOutCol, -X[pivotOutCol])
//    }
//
//
////    fun pivotAt(pivots: List<PivotPoint>) {
////        pivots.forEach { pivot -> applyPivot(pivot) }
////
////        val pivotColsByRow = IntArray(M.nRows) { -1 }
////        pivots.forEach { pivotColsByRow[it.row] = it.col }
////        for(i in pivotColsByRow.indices) {
////            val j = pivotColsByRow[i]
////            if(j != -1) {
////                if(i < pivotPoints.size) { // replacing an existing pivot
////                    val pivotedOutCol = pivotPoints[i]
////                    pivotPoints[i] = j
////                    semiPivot(pivotedOutCol, -X[pivotedOutCol])
////                } else { // add newly pivoted row to the block for pivoted rows
////                    if (i != pivotPoints.size) swapRows(i, pivotPoints.size)
////                    pivotPoints.add(j)
////                    X[j] = Fraction(B[i],M[i,j])
////                }
////            }
////        }
////    }
//
//
//    // Does the actual pivoting of M
//    // For a pivot at point (i,j), the k'th row of M and B is updated according to
//    // M_k' = (M_ij*M_k - M_kj*M_i)/G
//    // where G is the greatest common divisor of M_ij and M_kj
//    private fun applyPivot(pivot: PivotPoint) {
//        var Mij = M[pivot.row,pivot.col]
//        if (Mij < 0) {
//            M.timesAssignRow(pivot.row, -1)
//            Mij = -Mij
//        }
//
//        val col = M.columns[pivot.col].filter { it.key != pivot.row } // separate to prevent concurrent modification
//        for ((k, Mkj) in col) {
//            val gcd = gcd(Mij, Mkj)
//            M.timesAssignRow(k, Mij/gcd)
//            M.weightedRowPlusAssign(k, pivot.row, -Mkj/gcd)
//        }
//        basicColsByRow[pivot.row] = pivot.col
//    }
//
//
//    // perturbs the value of the unpivoted X elements, and updates the value
//    // of the pivoted X elements accordingly
////    fun columnPerturb() {
////        val perturbationDistribution = findAllPositivePivotCols().map {
////            Pair(it, 1.0)//exp(-changeInObservationNorm(it.value).toDouble()))
////        }
////        println("Found ${perturbationDistribution.size} possible perturbations")
////        val chosenCol = EnumeratedDistribution(perturbationDistribution).sample()
////        // update X_pivot for new X_unpivot value
////        applyNullVector(chosenCol, 1)
////    }
//
//
//    // perturbs X at the j'th element (col of M) by the given perturbation
//    // and updates the value of the basic X elements accordingly
//    // If p is the perturbation then the perturbation in X is
//    // X_p(i)' = -p*M_ij/M_ip
//    // X_j' = p
//    // where M_ip is the value of the currently pivoted-in element on the i'th row
//    fun applyNullVector(j: Int, perturbation: Fraction = maxPerturbation(j)) {
//        assert(!isBasicColumn(j))
//        X[j] += perturbation
//        M.columns[j].forEach { (i, Mij) ->
//            val pi = basicColsByRow[i]
//            X[pi] -= (perturbation * Mij)/M[i,pi]
//        }
//    }
//
//    fun applyNullVector(j: Int, perturbation: Int) {
//        applyNullVector(j, Fraction(perturbation,1))
//    }
//
//
////    fun pivotPerturb() {
////        val possiblePivots = findAllValidPivotPoints().map {
//////            Pair(it, 10.0.pow(-changeInObservationNorm(M[it.col]).toDouble()))
////            Pair(it, if(changeInObservationNorm(M.columns[it.col]) == 0) 1.0 else 0.0)
////        }
////        println("Found ${possiblePivots.size} possible pivot point. Prob sum is ${possiblePivots.map {it.second}.sum()}")
////        val choice = EnumeratedDistribution(possiblePivots).sample()
////        pivotAt(choice)
////        // setPivotedSolutionValues()
////    }
//
//
////    fun hybridPerturb() {
////        val xPerturbations = positivePerturbations()
////        val perturbationDistribution = (0 until M.nCols).asSequence()
////            .flatMap { j ->
////                xPerturbations[j]
////                    ?.asSequence()
////                    ?.mapNotNull { dXj -> if(dXj != 0) Pair(AbstractMap.SimpleEntry(j, dXj), 1.0) else null }
////                    ?:emptySequence()
////            }
////            .toList()
////        println("Found ${perturbationDistribution.size} possible perturbations")
////        val chosenPerturbation = EnumeratedDistribution(perturbationDistribution).sample()
////        perturbCol(chosenPerturbation.key, chosenPerturbation.value)
////    }
//
//
////    fun mcmcPerturb() {
////        val perturbationDistribution = this.perturbationDistribution?:mcmcTransitionProbs()
////        println("Found ${perturbationDistribution.size} possible perturbations")
////        val chosenPerturbation = perturbationDistribution.sample()
////        perturbCol(chosenPerturbation.column, chosenPerturbation.dx)
////        // choose whether to accept perturbation
////        val newPerturbationDistribution = mcmcTransitionProbs()
////        val acceptance =  perturbationDistribution.sumOfProbabilities() / newPerturbationDistribution.sumOfProbabilities()
////        if(acceptance < Random.nextDouble()) { // reject
////            perturbCol(chosenPerturbation.column, -chosenPerturbation.dx)
////            this.perturbationDistribution = perturbationDistribution
////        } else {
////            this.perturbationDistribution = newPerturbationDistribution
////        }
////    }
//
//
//    // Minimise CX
//    // Subject to:
//    // MX = B
//    // and
//    // X >= 0 and X is integer
//    // using Google OR-Tools
////    fun MPsolve(C: List<Double> = DoubleArray(M.nCols) {0.0}.asList()) {
////        M.IPsolve(B, C, "==").forEachIndexed { i, xi ->
////            X[i] = xi
////        }
////    }
//
//
//    class Perturbation(val column: Int, val dx: Int, val probability: Double, val cumulativeProb: Double)
//    // returns the probability of perturbations of each unpivoted column,
//    // given the current values of all other semi-pivoted columns
////    fun mcmcTransitionProbs(): List<Perturbation> {
////        var totalProb = 0.0
////        val perturbations = (0 until M.nCols).asSequence()
////            .flatMap { j ->
////                if (isBasicColumn(j)) emptySequence() else {
////                    val colProb = columnProbability(j)
////                    M.columns[j]
////                        .fold(IntRange(-X[j].toInt(), Int.MAX_VALUE)) { range, (i, Mij) ->
////                            val limit = X[pivotPoints[i]].toInt() / Mij
////                            if (Mij > 0 && limit < range.last) IntRange(range.start, limit)
////                            else if (Mij < 0 && limit > range.first) IntRange(limit, range.last)
////                            else range
////                        }
////                        .asSequence()
////                        .mapNotNull { dx ->
////                            if(dx != 0) {
////                                val prob = colProb.pow(0.5*dx) // square root to ensure acceptance equals ratio of sums
////                                totalProb += prob
////                                Perturbation(j, dx, prob, totalProb)
////                            } else null
////                        }
////                }
////            }.toList()
//////        perturbations.forEach { perturbation -> // normalise probabilities
//////            perturbation.probability /= totalProb
//////            perturbation.cumulativeProb /= totalProb
//////        }
////        return perturbations
////    }
//
//    fun List<Perturbation>.sample(): Perturbation {
//        val targetProb = Random.nextDouble(0.0, this.sumOfProbabilities())
//        val index = binarySearch(0, size) { (it.cumulativeProb - targetProb).sign.toInt() }
//        return this[index.absoluteValue]
//    }
//
//    fun List<Perturbation>.sumOfProbabilities(): Double = this.lastOrNull()?.cumulativeProb?:0.0
//
//        // returns the probability of this column
//    // i.e. the product of its event probabilities raised to the power of their occupation number
//    fun columnProbability(j: Int): Double {
//        return M.columns[j].fold(1.0) { prob, entry ->
//            prob * (eventProbabilities[entry.key]?.pow(entry.value)?:0.0)
//        }
//    }
//
//
//    // Sets column j to the given value,
//    // doing a full pivot where possible or semi-pivot otherwise
////    fun perturbCol(j: Int, perturbation: Int) {
////        assert(!isBasicColumn(j))
////        M.columns[j].find { (i, Mij) ->  // try to find a possible pivot point
////            Mij.absoluteValue == 1 && i < pivotPoints.size && X[pivotPoints[i]] == Fraction(perturbation,1) // TODO: Not convinced this is sensible
////        }?.run {
////            pivotAt(PivotPoint(key, j))
////        }?:run {
////            applyNullVector(j, perturbation)
////        }
////    }
//
//
//    // all state-changing pivot points that will not produce negative act values
////    fun findAllValidPivotPoints(): List<PivotPoint> {
////        return B.flatMap { Bentry ->
////            if(Bentry.key < pivotPoints.size && Bentry.value > 0) {
////                M.rows[Bentry.key].mapNotNull { MTentry ->
////                    val pivot = PivotPoint(Bentry.key, MTentry.key)
////                    if (isValidPivotPoint(pivot)) pivot else null
////                }
////            } else emptyList()
////        }
////    }
//
//
////    // finds all pairs of columns that when added together
////    //
////    fun findAllPerturbationColumnPairs() {
////
////    }
//
////    fun removeNegatives() {
////            while(B.entries.find { it.value < 0 }?.apply {
////                    val perturbations = findPositivePivotCols(key)
////                        .sortedBy { changeInObservationNorm(it.value) }
////                    val chosenCol = perturbations[0]
////                    X[chosenCol.index] += 1
////                    B -= chosenCol.value
////                    setPivotedSolutionValues(false)
////                    println("added col ${chosenCol.value}")
////                    println("observation norm is ${observationNorm()}")
////                    println("new solution is $X")
////            } != null) {}
////    }
//
//
////    fun findInitialSolution() {
////        val initialSolution = M.IPsolve(B, constraintType = "=")
////        val Xinit = SparseColIntMatrix.SparseIntColumn(initialSolution)
////        assert(M*Xinit == B)
////        Xinit.keys.forEach { M[it].keys.forEach { assert(it < pivotPoints.size) } }
////        for(j in initialSolution.indices) {
////            if(initialSolution[j] != 0) {
////                if(M[j].sparseSize > 1) {
////                    println("Pivoting...")
////                    val pivot = PivotPoint(
////                        M[j].entries.find { it.value == 1 && it.key < pivotPoints.size }?.key?:throw(RuntimeException("Can't pivot solution. Odd!"))
////                        ,j
////                    )
////                    applyPivot(pivot)
////                    pivotPoints[pivot.row] = pivot.col
////                    M = MT.transpose()
////                    println("Observation error ${observationNorm()}")
////                }
////            }
////        }
////        setPivotedSolutionValues()
////        // TODO: just a sanity check:
////        assert(M*X == B)
////    }
//
//
//    // pivots in all non-zero columns in the supplied solution
//    // then greedily pivots in any remaining rows
//    fun pivotToSolution(solution: SparseIntVector) {
//        for(entry in solution) {
//            if(!isBasicColumn(entry.key)) {
//                val pivot = PivotPoint(
//                    M.columns[entry.key]
//                        .filter { basicColsByRow[it.key] == -1 }
//                        .minBy { it.value.absoluteValue }
//                        ?.key
//                        ?:throw(RuntimeException("Can't pivot to solution. Must be a non-extreme solution (i.e. contains loops)!"))
//                    ,entry.key
//                )
//                applyPivot(pivot)
//            }
//        }
//        greedyPivotAll()
//    }
//
//
//    // pivots in any rows that aren't already pivoted in
//    // using a greedy algorithm:
//    // for each row in turn, choose to pivot on the column with smallest
//    // sparse size among all that have a +-1 in this row
//    fun greedyPivotAll() {
//        for(pivotRow in 0 until M.nRows) {
//            if(basicColsByRow[pivotRow] == -1) {
//                val pivotCol = M.rows[pivotRow]
//                    .minBy { it.value.absoluteValue*1024 + M.columns[it.key].sparseSize }
//                    ?.key
//                    ?: throw(RuntimeException("No elements to pivot on, equations must be degenerate"))
//                applyPivot(PivotPoint(pivotRow, pivotCol))
//            }
//        }
//    }
//
//
////    fun pivotOutNegatives() {
////        var tolerance = 0
////        while(negativeActCount() > 0) {
////            println("negative score is ${negativeActCount()}")
////            allPivotPoints()
////                .find { changeInNegativeActCountOnPivot(it) < tolerance }
////                ?.run {
////                    println("Pivoting out negative at $this")
////                    val expectedChange = changeInNegativeActCountOnPivot(this)
////                    println("expected count change is $expectedChange")
////                    val startScore = negativeActCount()
//////                    assert(M*X == B)
////                    pivotAt(this)
////                    println("actual change is ${negativeActCount() - startScore}")
//////                    assert(M*X == B)
////                    assert(negativeActCount() - startScore == expectedChange )
////                    tolerance = 0
////                }
////                ?:run {
////                    println("increasing tolerance")
////                    tolerance++
////                    if(tolerance == 3) throw(RuntimeException("No possible pivots to pivot out negatives"))
////                }
////        }
////        // setPivotedSolutionValues()
////    }
//
////    fun negativeActCount() = X.values.sumBy { if(it < Fraction.ZERO) -(it.toInt()) else 0 }
////        //B.sumBy { if(it.key < pivotPoints.size && it.value < 0) 1 else 0  }
//
//
//    // calculates the change in negativeActCount upon pivot at the given point
////    fun changeInNegativeActCountOnPivot(pivot: PivotPoint): Int {
////        val Bpiv = B[pivot.row]
////    //    assert(Bpiv >= 0)
////        val Mpiv = M[pivot.row, pivot.col]
////        assert(Mpiv.absoluteValue == 1)
////        return M.columns[pivot.col].sumBy { entry ->
////            if(entry.key != pivot.row && entry.key < pivotPoints.size) {
////                val bi = B[entry.key]
////                max(0, Bpiv*entry.value/Mpiv - bi) - max(0, -bi)
////            } else 0
////        } + if(Mpiv < 0) Bpiv else 0
////    }
//
//
//
//
//
//
//    // returns all points that have absolute value one and B value not equal to zero
////    fun allPivotPoints(): Sequence<PivotPoint> {
////        return (0 until pivotPoints.size).asSequence()
////            .flatMap { row ->
////                if(B[row] != 0) {
////                    M.rows[row].asSequence().mapNotNull { (col, Mij) ->
////                        if(Mij.absoluteValue == 1 && col != pivotPoints[row]) PivotPoint(row, col) else null
////                    }
////                } else emptySequence()
////            }
////    }
//
//    // finds pivot cols that are non-zero in given row
////    fun findPositivePivotCols(row: Int) =
////        M.columns
////            .withIndex()
////            .filter {
////                it.value.sparseSize > 1 && it.value.keys.contains(row) && isPositiveColumn(it.value) //&& changeInObservationNorm(it.value) <= 0
////            }
//
//
//    // returns the range of possible perturbations of each unpivoted column,
//    // given the current values of all other semi-pivoted columns
//    // or null if the column is already pivoted
////    fun positivePerturbations(): List<IntRange?> {
////        return Array(M.nCols) { j ->
////            if(isBasicColumn(j)) null else {
////                M.columns[j].asSequence()
////                    .fold(IntRange(-X[j].toInt(), Int.MAX_VALUE)) { range, (i, Mij) ->
////                        val limit = X[pivotPoints[i]].toInt()/Mij
////                        if(Mij > 0 && limit < range.last) IntRange(range.start, limit)
////                        else if(Mij < 0 && limit > range.first) IntRange(limit, range.last)
////                        else range
////                    }
////            }
////        }.asList()
////    }
//
//
//
//
//    // true if this column is currently pivoted
//    fun isBasicColumn(j: Int): Boolean {
//        val col = M.columns[j]
//        if(col.sparseSize != 1) return false
//        return basicColsByRow[col.first().key] == j
//    }
//
//
////    fun findAllPositivePivotCols(): List<Int> =
////        (0 until M.nCols).filter {
////            val col = M.columns[it]
////            col.sparseSize > 1 && isPositiveColumn(col)
////        }
//
////        M.columns.asSequence()
////            .withIndex()
////            .filter {
////                it.value.sparseSize > 1 && isPositiveColumn(it.value) //&& changeInObservationNorm(it.value) <= 0
////            }
////            .toList()
//
//
//
//
//    fun score(x: SparseIntVector): Double {
//        val rootOffset = M.rows[M.nRows-1].dotProd(x).absoluteValue
//        return -(x.count {it.key != 1-M.nRows && (it.key < 0 || it.value < 0) } + rootOffset).toDouble()
//    }
//
//    // true iff subtracting this column from B increases the 1-norm of the observation
////    fun changeInObservationNorm(col: SparseIntVector): Int {
////        var dErr = 0
////        for(entry in col) {
////            if(entry.key >= pivotPoints.size) {
////                val bi = B[entry.key]
////                dErr += (bi - entry.value).absoluteValue - bi.absoluteValue
////            }
////        }
////        return dErr
////    }
//
//    // returns the 1-norm of the observation
////    fun observationNorm(): Int {
////        var norm = 0
////        for(entry in B) {
////            if(entry.key >= pivotPoints.size) {
////                norm += entry.value.absoluteValue
////            }
////        }
////        return norm
////    }
//
//
//
//    // true iff setting this column to one will not create any negative values in B
////    fun isPositiveColumn(col: SparseIntVector): Boolean {
////        return col.all { entry ->
//////                B[entry.key] >= entry.value
////                entry.key < pivotPoints.size && X[pivotPoints[entry.key]].numerator >= entry.value
////            }
////    }
//
//    // true iff this point equals +-1, is not already pivoted-in and
//    // pivoting at this point will not create any negative values in B
////    fun isValidPivotPoint(pivot: PivotPoint): Boolean {
////        val pivCol = M.columns[pivot.col]
////        if(pivCol[pivot.row].absoluteValue != 1 || pivCol.sparseSize <= 1) return false
////        val k = B[pivot.row]/pivCol[pivot.row]
////        return k >= 0 && pivCol.all { B[it.key] >= k*it.value || it.key >= pivotPoints.size }
////    }
//
//
////    fun swapRows(i1: Int, i2: Int) {
////        if((i1 < pivotPoints.size) xor (i2 < pivotPoints.size)) throw(RuntimeException("Can't swap pivoted row with non-pivoted row"))
////        if(i1 == i2) return
////        M.swapRows(i1,i2)
////        val tmp = B[i1]
////        B[i1] = B[i2]
////        B[i2] = tmp
////        if(i1 < pivotPoints.size) {
////            val ppTmp = pivotPoints[i1]
////            pivotPoints[i1] = pivotPoints[i2]
////            pivotPoints[i2] = ppTmp
////        }
////    }
//
//    fun calcB(): SparseIntVector {
//        val B = HashIntVector()
//        for(i in 0 until M.nRows) {
//            B[i] = X.sumByDouble { (j, Xj) -> (Xj*M[i,j]).toDouble() }.roundToInt()
//        }
//        return B
//    }
//
//    // returns M with B added as a rightmost column
//    fun toCompositeMatrix(): HashRowColIntMatrix {
//        val compositeMatrix = M.copy()
//        compositeMatrix.resize(M.rows.size, M.columns.size + 1)
//        compositeMatrix.setColumn(M.columns.size, calcB())
//        return compositeMatrix
//    }
//
//    override fun toString(): String {
//        val string = StringBuilder()
//        val activeCols = CharArray(M.nCols*3) { '.' }
//        basicColsByRow.forEach { j ->
//            activeCols[j*3 + 1] = '#'
//        }
//        string.append(activeCols)
//        string.append('\n')
//        string.append(toCompositeMatrix().toString())
//        return string.toString()
//    }
//
//    fun toSparsityString(): String {
//        val string = StringBuilder()
//        val activeCols = CharArray(M.nCols) { '.' }
//        basicColsByRow.forEach { j ->
//            activeCols[j] = '#'
//        }
//        string.append(activeCols)
//        string.append('\n')
//        string.append(toCompositeMatrix().toSparsityString())
//        return string.toString()
//
//    }
//
//
//    // Adds a 'fake root node that has value zero and edges from it to each
//    // source node, so that the whole matrix becomes a single tree with only forward edges
//    // takes a list of row indices which identify the source rows
//    // and adds the fake root to the bottom of this matrix (increasing nRows by 1)
//    // and a new action for each source (increasing nCols by the number of sources)
//    // returns the row index of the added root
////    fun addSourceRoot(sources: List<Int>): Int {
////        val i = M.nRows
////        val j = M.nCols
////        M.resize(i + 1, j + sources.size)
////        for(col in sources.indices) {
////            M[i, j + col] = -1
////            M[sources[col], j + col] = 1
////        }
////        return M.nRows-1
////    }
//
//
//    data class PivotPoint(val row: Int, val col: Int)
//
////    enum class RowType {
////        PIVOTED,
////        PIVOTEDROOT,
////        UNPIVOTED,
////        UNPIVOTEDROOT
////
////        fun isPivoted() = (this == PIVOTED || this == PIVOTEDROOT)
////    }
//
////    inner class TreeConstructor {
////        val ROOT = -2
////        val distanceUpperBounds =
////            HashMap<Int, MutableList<PivotPoint>>() // upper bounds of distance to root. Maps distance to pivot points
////        val colDistances: IntArray              // lower bound on distance when pivoting on this column
////                                                // when the number of undecided rows for a column becomes 1, this bound is tight
////        var distanceLowerBound = 0
////        var highestBound = 0
////
////        constructor(rootRows: List<Int>, initialPivots: Map<Int,Int> = emptyMap()) {
////            distanceUpperBounds[0] = rootRows.map { PivotPoint(it, ROOT) }.toMutableList()
////
////            val Mcols = HashColIntMatrix(M)
////            colDistances = IntArray(M.nCols) { 1 }
////            while(distanceLowerBound <= highestBound) {
////                val leaves = distanceUpperBounds[distanceLowerBound]
////                    ?.filter { it.col == ROOT || Mcols[it.row, it.col] != 0 } // remove if already reduced
////                    ?.distinctBy { it.row }
////                    ?:listOf()
////                distanceUpperBounds[distanceLowerBound] = leaves.toMutableList()
////                leaves
////                    .forEach { pivotPoint -> // recalculate colDistances and add new children to upper bounds
////                        M.rows[pivotPoint.row]
////                            .forEach { (childCol, childVal) ->
////                                assert(Mcols[pivotPoint.row, childCol] != 0)
////                                val priorSize = Mcols.columns[childCol].sparseSize
//////                                println("initial val = ${M[pivotPoint.row, childCol]}")
////                                Mcols[pivotPoint.row, childCol] = 0 // use Mcols to keep track of followed values
//////                                println("prior size = $priorSize new size = ${M.columns[childCol].sparseSize}")
////                                assert(Mcols.columns[childCol].sparseSize == priorSize - 1)
////                                colDistances[childCol] += distanceLowerBound
////                                if(Mcols.columns[childCol].sparseSize == 1 && childVal < 0) { // follow children forward in time when reduced
////                                    distanceUpperBounds
////                                        .getOrPut(colDistances[childCol]) {ArrayList()}
////                                        .add(PivotPoint(Mcols.columns[childCol].keys.first(), childCol))
////                                    if(colDistances[childCol] > highestBound) highestBound = colDistances[childCol]
////                                }
////                            }
////                    }
////                distanceLowerBound++
////            }
////
//////            assert(highestBound <= distanceLowerBound)
////
//////            println("Pivot map is $distanceUpperBounds")
////            // now we can read off the tree from the distanceUpperBounds
////            // pivoting in reverse order is most efficient
////
////            val pivotsInOrder = (distanceLowerBound-1 downTo 1)
////                .flatMap {
////                    distanceUpperBounds[it]
////                        ?.map { originalPivot ->
////                            initialPivots[originalPivot.row]
////                                ?.run { PivotPoint(originalPivot.row, this) }
////                                ?:originalPivot
////                        }
////                        ?:emptyList()
////                }
////
////            pivotAt(pivotsInOrder)
////            greedyPivotAll() // finally, pivot on tree roots
////        }
////
////    }
//}