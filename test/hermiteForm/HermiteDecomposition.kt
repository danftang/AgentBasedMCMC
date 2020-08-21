package hermiteForm

class HermiteDecomposition {
    val mask: BooleanArray // true if act with this index is the head of a basis
//    val pseudoInverse: SparseColIntMatrix
    val hermiteForm: SparseColIntMatrix
    val U: SparseColIntMatrix
    val nullBasis: SparseColIntMatrix
    val baseState: SparseColIntMatrix.SparseIntColumn

    constructor(abm: SparseColIntMatrix, B: SparseColIntMatrix.SparseIntColumn) {
        hermiteForm = SparseColIntMatrix(abm)
        U = hermiteForm.hermiteDecomposition()
        nullBasis = SparseColIntMatrix(U.subList(abm.nRows, U.nCols), U.nRows)
        nullBasis.upperTriangularise()
        nullBasis.upperTriangularThin()

        mask = BooleanArray(nullBasis.nRows) { false }
        for(col in nullBasis.columns) {
            val basisAct = col.keys.max()
            if(basisAct != null) mask[basisAct] = true
        }

        baseState = baseSolution(B)
    }

    // Transforms an act vector by removing all elements that correspond to elements in the basis vector
    fun actToNonBasisVector(actVector: SparseColIntMatrix.SparseIntColumn): SparseColIntMatrix.SparseIntColumn {
        var nonBasisIndex = 0
        val nonBasisVector = SparseColIntMatrix.SparseIntColumn()
        for(i in mask.indices) {
            if(!mask[i]) {
                val actVal = actVector[i]
                if(actVal != 0) nonBasisVector[nonBasisIndex] = actVal
                nonBasisIndex++
            }
        }
        return nonBasisVector
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
        val baseHermiteVector = hermiteForm.solveLowerTriangular(B)
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
    fun toConstraintMatrix(): SparseColIntMatrix {
        val constraints = nullBasis.transpose()
        constraints.removeAll(constraints.filterIndexed { i,_ -> mask[i]})
        return constraints.transpose()
    }

}