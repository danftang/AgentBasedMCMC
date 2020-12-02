package lib.sparseIntMatrix

import kotlin.math.min

interface SparseColIntMatrix: SparseIntMatrix {
    val columns:    List<SparseIntVector>

    override val entries: Sequence<SparseIntMatrix.Entry>
            get() = columns
                .asSequence()
                .mapIndexed { j, col ->
                    col.asSequence().map { SparseIntMatrix.Entry(it.key, j, it.value) }
                }
                .flatten()

    fun clearColumn(j: Int)

    fun swapCols(j1: Int, j2: Int)

    fun setColumn(col: Int, vector: Iterable<Map.Entry<Int,Int>>) {
        clearColumn(col)
        for(entry in vector) {
            this[entry.key, col] = entry.value
        }
    }

    fun weightedColPlusAssign(i1: Int, i2: Int, weight: Int) {
        for(entry in columns[i2]) {
            plusAssign(entry.key, i1, weight*entry.value)
        }
    }


    override operator fun times(X: SparseIntVector): HashIntVector {
        val result = HashIntVector()
        for(entry in X) {
            result.weightedPlusAssign(columns[entry.key], entry.value)
        }
        return result
    }


    override operator fun times(X: IntArrayVector): IntArrayVector {
        val result = IntArrayVector(nRows)
        for(i in X.indices) {
            for(entry in columns[i]) {
                result[entry.key] += entry.value * X[i]
            }
        }
        return result
    }


    fun removeColumns(colsToRemove: Iterable<Int>)

    fun replaceNonZeroElementsInCol(col: Int, map: (row: Int, value: Int) -> Int)


    fun hermiteDecomposition(): HashColIntMatrix {
        val U = HashColIntMatrix.identity(this.nCols)
        for(eqnIndex in 0 until this.nRows) {
            var indexColMagnitude = this[eqnIndex,eqnIndex]
            for(colIndex in this.nCols-1 downTo eqnIndex+1) {
                val swapColMagnitude = this[eqnIndex,colIndex]
                if(swapColMagnitude != 0) {
                    if(indexColMagnitude != 0) {
                        val g = HashColIntMatrix.ExtendedEuclid(indexColMagnitude, swapColMagnitude)
                        val indexMultiplier = swapColMagnitude/g.gcd
                        val swapColMultiplier = indexColMagnitude/g.gcd
                        hermiteColumnSwap(eqnIndex, colIndex, g.y, swapColMultiplier, indexMultiplier)
                        hermiteColumnSwap(eqnIndex, colIndex, g.y, swapColMultiplier, indexMultiplier)
                        indexColMagnitude = this[eqnIndex,eqnIndex]
                    } else {
                        this.swapCols(eqnIndex, colIndex)
                        U.swapCols(eqnIndex, colIndex)
                        indexColMagnitude = swapColMagnitude
                    }
                }
            }
            // now reduce to the left
//            for(colIndex in 0 until eqnIndex) {
//                val swapColMagnitude = this[eqnIndex,colIndex]
//                if(swapColMagnitude != 0) {
//                    val reductionWeight = -swapColMagnitude / indexColMagnitude
//                    this[colIndex].weightedPlusAssign(indexCol, reductionWeight)
//                    U[colIndex].weightedPlusAssign(UindexCol, reductionWeight)
//                }
//            }
        }
        return U
    }


    // sets
    // c1' = X*c1 + Y*c2
    // c2' = B*c1 - A*c2
    // where X is an integer, A != 0 and XA + YB = 1
    // which is equivalent to
    // c2' = B*c1 - A*c2
    // c1' = (c1 - Y*c2')/A
    private fun hermiteColumnSwap(col1: Int, col2: Int, Y: Int, A: Int, B: Int) {
        replaceNonZeroElementsInCol(col2) {_,v -> -A*v}
        weightedColPlusAssign(col2, col1, B)
        weightedColPlusAssign(col1, col2, -Y)
        replaceNonZeroElementsInCol(col1) {_,v -> v/A}
    }


    // solves LX=Y where L is the lower triangular
    // part of this matrix
    // Assumes that L^-1 is integer
    fun solveLowerTriangular(Y: SparseIntVector): HashIntVector {
        val X = HashIntVector(Y)
        for(i in 0 until min(nCols, nRows)) {
            val Xi = X[i]/this[i,i]
//            assert(Xi*this[i,i] == X[i])
            X[i] = Xi
            for(entry in columns[i]) {
                if(entry.key > i) X[entry.key] -= Xi*entry.value
            }
        }
        return X
    }

}