
import lib.SparseColIntMatrix
import java.util.*
import kotlin.math.absoluteValue

class HermiteDecomposition {
//    val mask: BooleanArray // true if act with this index is the head of a null basis
    val nullMaskMatrix: SparseColIntMatrix // pre-multiplication with this matrix removes rows in nullBasis that have a single 1 in them
//    val pseudoInverse: SparseColIntMatrix
    val H: SparseColIntMatrix
    val U: SparseColIntMatrix
    val nullBasis: SparseColIntMatrix
    val solutionBasis: SparseColIntMatrix
    val baseState: SparseColIntMatrix.SparseIntColumn

    constructor(abm: SparseColIntMatrix, observation: SparseColIntMatrix.SparseIntColumn) {
        H = SparseColIntMatrix(abm)
        U = H.hermiteDecomposition()
        nullBasis = U.subMartixView(abm.nRows, U.nCols)
        nullBasis.upperTriangularise()
        nullBasis.upperTriangularThin()
        solutionBasis = U.subMartixView(0, abm.nRows)

//        mask = BooleanArray(nullBasis.nRows) { false }
//        for(col in nullBasis.columns) {
//            val basisAct = col.keys.max()
//            if(basisAct != null) mask[basisAct] = true
//        }
        nullMaskMatrix = calculateMaskMatrix()

        baseState = baseSolution(observation)
    }


    // Tree decomposition
    constructor(abm: SparseColIntMatrix, observation: SparseColIntMatrix.SparseIntColumn, treeRootRow: Int) {
        val decomp = TreeConstructor(abm.copy(), treeRootRow)
        H = decomp.H
        U = decomp.U
        nullBasis = U.subMartixView(abm.nRows, U.nCols)
        solutionBasis = U.subMartixView(0, abm.nRows)
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


    private fun calculateMaskMatrix(): SparseColIntMatrix {
        val M = SparseColIntMatrix(nullBasis.nRows, U.nRows)
        var row = 0
        nullBasis
            .transpose()
            .mapIndexedNotNull { i, col ->
                if(col.sparseSize > 1) i else null
            }
            .forEach {i ->
                M[row++,i]=1
            }
        M.nRows = row
        return M
    }

    fun basisVectorToActs(basisVector: SparseColIntMatrix.SparseIntColumn): SparseColIntMatrix.SparseIntColumn {
        return baseState + nullBasis * basisVector
    }

    fun basisVectorToActs(basisVector: IntArray): SparseColIntMatrix.SparseIntColumn {
        return baseState + nullBasis * basisVector
    }

    fun basisVectorToActs(basisVector: List<Int>): SparseColIntMatrix.SparseIntColumn {
        return baseState + nullBasis * basisVector
    }


    // calculate the base solution, X, to the original abm constraints AX = B
    // The base solution is the one that maps to the zero null basis vector
    fun baseSolution(B: SparseColIntMatrix.SparseIntColumn): SparseColIntMatrix.SparseIntColumn {
        val baseHermiteVector = H.solveLowerTriangular(B)
        val baseState = U * baseHermiteVector
        // now remove surplus null vectors so all null vectors so negative null basis additions are discounted
        for(basis in nullBasis.asReversed()) {
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
        val HT: SparseColIntMatrix
        val H: SparseColIntMatrix
        val U: SparseColIntMatrix
        var basisSize = 1
        val rootRow: Int
        val colIsInteraction: Array<Boolean>            // true if this action is an interaction

        constructor(H: SparseColIntMatrix, rootRow: Int) {
            this.H = H
            this.rootRow = rootRow
            U = SparseColIntMatrix.Identity(H.nCols)
            pivotPoints = ArrayDeque<Pair<Int,Int>>(H.nRows)    // (row,col)
            HT = H.transpose()
            colIsInteraction = Array(H.nCols) { j ->
                H[j].sparseSize > 2
            }
            HT[rootRow].keys
                .filter { j -> H[j].sparseSize == 2 }
                .forEach { j ->
                    addPivotPoint(H[j].keys.find {it != rootRow}!! ,j)
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
                val col = H[j]
                if(col.sparseSize == 1 && col.keys.first() == rootRow) {
                    if(U[j].values.sumBy { it.absoluteValue } < minReductionSize) {
                        minReductionSize = U[j].sparseSize
                        minReducerColumn = j
                    }
                }
            }
            if(minReducerColumn != -1) {
                println("reducing root on column $minReducerColumn ${H[minReducerColumn]}")
                val indexMagnitude = H[rootRow,minReducerColumn]
                val pivotCol = H[minReducerColumn]
                val UPivotCol = U[minReducerColumn]
                assert(indexMagnitude.absoluteValue == 1)
                for(j in 0 until H.nCols) {
                    if(j != minReducerColumn) {
                        val swapMagnitude = H[rootRow, j]
                        if (swapMagnitude != 0) {
                            val weight = -swapMagnitude / indexMagnitude
                            H[j].weightedPlusAssign(pivotCol, weight)
                            U[j].weightedPlusAssign(UPivotCol, weight)
                        }
                    }
                }
            }

            // now rearrange columns into identity
            for(j in 0 until H.nCols) {
                var swapj = H[j].keys.firstOrNull()
                while(swapj != null && swapj != j) {
                    H.swapCols(j, swapj)
                    U.swapCols(j, swapj)
                    swapj = H[j].keys.firstOrNull()
                }
            }

            // now ensure all columns are +ve
            for(j in 0 until H.nRows) {
                if(H[j,j] < 0) {
                    H[j] *= -1
                    U[j] *= -1
                }
            }
        }

        fun reduceOn(row: Int, col: Int) {
            val pivotVal = H[row,col]
//            println("reducing on $row $col ${H[col]}")
            assert(pivotVal.absoluteValue == 1)
            val HpivotCol = H[col]
            val UpivotCol = U[col]
            HT[row].entries.forEach { rowElement ->
                val colToReduce = rowElement.key
                if(colToReduce != col) {
                    val weight = -rowElement.value/pivotVal
                    H[colToReduce].weightedPlusAssign(HpivotCol, weight)
                    U[colToReduce].weightedPlusAssign(UpivotCol, weight)
                    val reducedSize = H[colToReduce].keys.size - if(H[colToReduce].keys.contains(rootRow)) 1 else 0
                    if(reducedSize == 1) {
                        if(colIsInteraction[colToReduce]) {
                            addPendingInteraction(colToReduce)
                        } else {
                            addPivotPoint(H[colToReduce].keys.find { it != rootRow }!!, colToReduce)
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
                .getOrPut(U[j].values.sumBy { it.absoluteValue } ) { ArrayList() }
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
            val pivotSize = U[pivotPoint.second].sparseSize
            if(pivotSize > basisSize) {
                basisSize = pivotSize
                pendingInteractions
                    .get(basisSize)
                    ?.forEach { j ->
                        val i = H[j].keys.find { it != rootRow }
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