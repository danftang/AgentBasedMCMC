
import lib.sparseIntMatrix.*
import java.util.*
import kotlin.math.absoluteValue

class HermiteDecomposition {
//    val mask: BooleanArray // true if act with this index is the head of a null basis
    val nullMaskMatrix: HashColIntMatrix // pre-multiplication with this matrix removes rows in nullBasis that have a single 1 in them
//    val pseudoInverse: SparseColIntMatrix
    val H: HashRowColIntMatrix
    val U: HashColIntMatrix
    val nullBasis: HashColIntMatrix
    val solutionBasis: HashColIntMatrix
    val baseState: HashIntVector

    constructor(abm: HashColIntMatrix, observation: HashIntVector) {
        H = HashRowColIntMatrix(abm)
        U = H.hermiteDecomposition()
        nullBasis = U.subMatrixView(abm.nRows, U.nCols)
        nullBasis.upperTriangularise()
        nullBasis.upperTriangularThin()
        solutionBasis = U.subMatrixView(0, abm.nRows)

//        mask = BooleanArray(nullBasis.nRows) { false }
//        for(col in nullBasis.columns) {
//            val basisAct = col.keys.max()
//            if(basisAct != null) mask[basisAct] = true
//        }
        nullMaskMatrix = calculateMaskMatrix()

        baseState = baseSolution(observation)
    }


    // Tree decomposition
    constructor(abm: HashColIntMatrix, observation: SparseIntVector, treeRootRow: Int) {
        val decomp = TreeConstructor(HashRowColIntMatrix(abm), treeRootRow)
        H = decomp.H
        U = decomp.U
        nullBasis = U.subMatrixView(abm.nRows, U.nCols)
        solutionBasis = U.subMatrixView(0, abm.nRows)
        baseState = solutionBasis * observation

        nullMaskMatrix = calculateMaskMatrix()
    }

    // Transforms an act vector by removing all elements that correspond to elements in the basis vector
//    fun actToNonBasisVector(actVector: SparseColIntMatrix.SparseIntColumn): SparseColIntMatrix.SparseIntColumn {
//        var nonBasisIndex = 0
//        val nonBasisVector = SparseColIntMatrix.SparseIntColumn()
//        for(i in mask.indices) {
//            if(!mask[i]) {
//                val actVal = actVector[i]
//                if(actVal != 0) nonBasisVector[nonBasisIndex] = actVal
//                nonBasisIndex++
//            }
//        }
//        return nonBasisVector
//    }


    private fun calculateMaskMatrix(): HashColIntMatrix {
        val M = HashColIntMatrix(nullBasis.nRows, U.nRows)
        var row = 0
        nullBasis
            .transpose()
            .columns.mapIndexedNotNull { i, col ->
                if(col.sparseSize > 1) i else null
            }
            .forEach {i ->
                M[row++,i]=1
            }
        M.nRows = row
        return M
    }

    fun basisVectorToActs(basisVector: HashIntVector): HashIntVector {
        return baseState + nullBasis * basisVector
    }

    fun basisVectorToActs(basisVector: IntArrayVector): HashIntVector {
        return basisVectorToActs(basisVector.toSparseIntVector())
    }

//    fun basisVectorToActs(basisVector: List<Int>): SparseColIntMatrix.SparseIntColumn {
//        return baseState + nullBasis * basisVector
//    }


    // calculate the base solution, X, to the original abm constraints AX = B
    // The base solution is the one that maps to the zero null basis vector
    fun baseSolution(B: HashIntVector): HashIntVector {
        val baseHermiteVector = H.solveLowerTriangular(B)
        val baseState = U * baseHermiteVector
        // now remove surplus null vectors so all null vectors so negative null basis additions are discounted
        for(basis in nullBasis.columns.asReversed()) {
            val maxIndex = basis.keys.max()
            if(maxIndex != null) {
                baseState.weightedPlusAssign(basis, -baseState[maxIndex]/basis[maxIndex])
            }
        }
        return baseState
    }

    // the null basis vector with the singleton basis rows removed
//    fun toConstraintMatrix(): SparseColIntMatrix {
//        val constraints = nullBasis.transpose()
//        constraints.removeAll(constraints.filterIndexed { i,_ -> mask[i]})
//        return constraints.transpose()
//    }


    // Constructs the decomposition by breadth first tree serach
    // starting with a given root row.
    //
    class TreeConstructor {
        val pendingInteractions = HashMap<Int,ArrayList<Int>>() // queue of interactions to possibly join tree. Maps basis size to column indices
        val pivotPoints: ArrayDeque<Pair<Int,Int>>      // Pivot points waiting to be reduced (row,col)
        val rowsInTree = HashSet<Int>()                 // rows already reduced
 //       val HT: HashColIntMatrix
        val H: HashRowColIntMatrix
        val U: HashColIntMatrix
        var basisSize = 1
        val rootRow: Int
        val colIsInteraction: Array<Boolean>            // true if this action is an interaction

        constructor(H: HashRowColIntMatrix, rootRow: Int) {
            this.H = H
            this.rootRow = rootRow
            U = HashColIntMatrix.identity(H.nCols)
            pivotPoints = ArrayDeque<Pair<Int,Int>>(H.nRows)    // (row,col)
            colIsInteraction = Array(H.nCols) { j ->
                H.columns[j].sparseSize > 2
            }
            H.rows[rootRow].keys
                .filter { j -> H.columns[j].sparseSize == 2 }
                .forEach { j ->
                    addPivotPoint(H.columns[j].keys.find {it != rootRow}!! ,j)
                } // add all non-interaction events
            while(!pivotPoints.isEmpty()) {
                val pivotPoint = getNextPivotPoint()
                reduceOn(pivotPoint.first, pivotPoint.second)
            }
            println("Reduced on ${rowsInTree.size} rows")

            // now reduce root
            var minReductionSize = Int.MAX_VALUE
            var minReducerColumn = -1
            for(j in 0 until H.nCols) {
                val col = H.columns[j]
                if(col.sparseSize == 1 && col.keys.first() == rootRow) {
                    if(U.columns[j].values.sumBy { it.absoluteValue } < minReductionSize) {
                        minReductionSize = U.columns[j].sparseSize
                        minReducerColumn = j
                    }
                }
            }
            if(minReducerColumn != -1) {
                println("reducing root on column $minReducerColumn ${H.columns[minReducerColumn]}")
                val indexMagnitude = H[rootRow,minReducerColumn]
                assert(indexMagnitude.absoluteValue == 1)
                for(j in 0 until H.nCols) {
                    if(j != minReducerColumn) {
                        val swapMagnitude = H[rootRow, j]
                        if (swapMagnitude != 0) {
                            val weight = -swapMagnitude / indexMagnitude
                            H.weightedColPlusAssign(j, minReducerColumn, weight)
                            U.weightedColPlusAssign(j, minReducerColumn, weight)
                        }
                    }
                }
            }

            // now rearrange columns into identity
            for(j in 0 until H.nCols) {
                var swapj = H.columns[j].keys.firstOrNull()
                while(swapj != null && swapj != j) {
                    H.swapCols(j, swapj)
                    U.swapCols(j, swapj)
                    swapj = H.columns[j].keys.firstOrNull()
                }
            }

            // now ensure all columns are +ve
            for(j in 0 until H.nRows) {
                if(H[j,j] < 0) {
                    H.replaceNonZeroElementsInCol(j) { _, v -> -v}
                    U.replaceNonZeroElementsInCol(j) { _, v -> -v}
                }
            }
        }

        fun reduceOn(row: Int, col: Int) {
            val pivotVal = H[row,col]
//            println("reducing on $row $col ${H[col]}")
            assert(pivotVal.absoluteValue == 1)
            H.rows[row].forEach { rowElement ->
                val colToReduce = rowElement.key
                if(colToReduce != col) {
                    val weight = -rowElement.value/pivotVal
                    H.weightedColPlusAssign(colToReduce, col, weight)
                    U.weightedColPlusAssign(colToReduce, col, weight)
                    val reducedSize = H.columns[colToReduce].sparseSize - if(H.columns[colToReduce].keys.contains(rootRow)) 1 else 0
                    if(reducedSize == 1) {
                        if(colIsInteraction[colToReduce]) {
                            addPendingInteraction(colToReduce)
                        } else {
                            addPivotPoint(H.columns[colToReduce].keys.find { it != rootRow }!!, colToReduce)
                        }
                    }
//                    val childMaxVal = H[colToReduce].values.maxBy { it.absoluteValue }?.absoluteValue?:-1
//                    if(reducedSize == 1 && childMaxVal == 1) {
//                        val newChildRow = H[colToReduce].keys.find { it != rootRow}!!
//                        if(!rowsInTree.contains(newChildRow)) {
//                            pivotPoints.addLast(Pair(newChildRow, colToReduce))
//                            rowsInTree.add(newChildRow)
//                        }
//                    }
                }
            }
        }


        fun addPendingInteraction(j: Int) {
            pendingInteractions
                .getOrPut(U.columns[j].values.sumBy { it.absoluteValue } ) { ArrayList() }
                .add(j)

        }

        fun addPivotPoint(i: Int, j: Int) {
            if(!rowsInTree.contains(i)) {
                pivotPoints.addLast(Pair(i, j))
                rowsInTree.add(i)
            }
        }

        fun getNextPivotPoint(): Pair<Int,Int> {
            val pivotPoint = pivotPoints.pollFirst()
            val pivotSize = U.columns[pivotPoint.second].sparseSize
            if(pivotSize > basisSize) {
                basisSize = pivotSize
                pendingInteractions
                    .get(basisSize)
                    ?.forEach { j ->
                        val i = H.columns[j].keys.find { it != rootRow }
                        if(i != null && !rowsInTree.contains(i)) {
                            pivotPoints.addLast(Pair(i, j))
                            rowsInTree.add(i)
                        }
                    }
                pendingInteractions.remove(basisSize)
            }
            return pivotPoint
        }

    }



    companion object {


//        fun unconditionalChildren(H: SparseColIntMatrix, parentRow: SparseColIntMatrix.SparseIntColumn, rootRow: Int) =
//            parentRow.keys
//                .filter { j ->
//                    val col = H[j]
//                    val reducedSize = col.keys.size - if(col.keys.contains(rootRow)) 1 else 0
//                    val childMaxVal = col.values.maxBy { it.absoluteValue }?.absoluteValue?:0
//                    reducedSize == 1 && childMaxVal == 1
//                }
//                .map { j ->
//                    Pair(H[j].keys.find { it != rootRow }!!, j)
//                }


    }
}