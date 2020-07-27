package hermiteForm

import org.apache.commons.math3.util.ArithmeticUtils

class SparseColIntMatrix: ArrayList<SparseColIntMatrix.SparseIntColumn> {
    val nRows: Int
    val nCols: Int
        get() = this.size

    constructor(nRows: Int, nCols: Int): super(nCols) {
        for(col in 1..nCols) add(SparseIntColumn())
        this.nRows = nRows
    }

    constructor(copy: SparseColIntMatrix): super(copy.nCols) {
        for(col in copy) add(SparseIntColumn(col))
        this.nRows = copy.nRows
    }

    operator fun set(row: Int, col: Int, element: Int) {
        this[col][row] = element
    }

    operator fun get(row: Int, col: Int): Int = this[col].getOrDefault(row,0)

    fun swapCols(col1: Int, col2: Int) {
        val tmp = this[col1]
        this[col1] = this[col2]
        this[col2] = tmp
    }

    fun setRow(row: Int, rowEntries: Collection<Map.Entry<Int,Int>>) {
        for(entry in rowEntries) {
            this[entry.key][row] = entry.value
        }
    }

    // ((Y+X)x -(Yx - Xy))/X = x + y
    // or ((Y+X)Yx -Y(Yx - Xy))/XY = x + y


    fun hermiteDecomposition(): SparseColIntMatrix {
        val U = Identity(this.nCols)
        for(eqnIndex in 0 until this.nRows) {
            var indexCol = this[eqnIndex]
            var UindexCol = U[eqnIndex]
            var indexColMagnitude = indexCol.getOrDefault(eqnIndex,0)
            for(colIndex in eqnIndex+1 until this.nCols) {
                val swapCol = this[colIndex]
                val UswapCol = U[colIndex]
                val swapColMagnitude = swapCol.getOrDefault(eqnIndex,0)
                if(swapColMagnitude != 0) {
                    if(indexColMagnitude != 0) {
                        val g = ExtendedEuclid(swapColMagnitude, indexColMagnitude)
                        val indexMultiplier = swapColMagnitude/g.gcd
                        val swapColMultiplier = indexColMagnitude/g.gcd
                        hermiteColumnSwap(indexCol, swapCol, g.y, swapColMultiplier, indexMultiplier)
                        hermiteColumnSwap(UindexCol, UswapCol, g.y, swapColMultiplier, indexMultiplier)
                        indexColMagnitude = indexCol.getOrDefault(eqnIndex,0)
                    } else {
                        this.swapCols(eqnIndex, colIndex)
                        indexCol = this[eqnIndex]
                        U.swapCols(eqnIndex, colIndex)
                        UindexCol = U[eqnIndex]
                        indexColMagnitude = swapColMagnitude
                    }
                }
            }
        }
        return U
    }

    // sets
    // c1' = X*c1 + Y*c2
    // c2' = B*c1 - A*c2
    // where X is an integer, A != 0 and XA + YB = 1
    // which is equivalent to
    // c2 = B*c1 - A*c2
    // c1 = (c1 - Y*c2)/A
    private fun hermiteColumnSwap(c1: SparseIntColumn, c2: SparseIntColumn, Y: Int, A: Int, B: Int) {
        c2 *= -A
        c2.weightedPlusAssign(c1,B)
//        c1.weightedPlusAssign(c2,-Y)
//        c1 /= A
    }



    override fun toString(): String {
        val out = StringBuilder()
        for(row in 0 until nRows) {
            for (col in this) {
                val v = col.getOrDefault(row,0)
                out.append(v)
                out.append(' ')
            }
            out.appendln()
        }
        return out.toString()
    }

    fun toSparsityString(): String {
        val out = StringBuilder()
        for(row in 0 until nRows) {
            for (col in this) {
                val v = col.getOrDefault(row,0)
                out.append(if(v == 0) '.' else '#')
            }
            out.appendln()
        }
        return out.toString()
    }


    // map from column index to non-zero entry value
    class SparseIntColumn: HashMap<Int,Int> {

        constructor(): super(2)

        constructor(copy: SparseIntColumn): super(copy)

        operator fun timesAssign(multiplier: Int) {
            for(entry in this.entries) {
                entry.setValue(entry.value * multiplier)
            }
        }

        operator fun divAssign(denominator: Int) {
            for(entry in this.entries) {
                entry.setValue(entry.value/denominator)
            }
        }

        operator fun minusAssign(otherCol: SparseIntColumn) {
            for(entry in otherCol.entries) {
                this.merge(entry.key, -entry.value) { a, b ->
                    val result = a + b
                    if(result != 0) result else null
                }
            }
        }

        // this += weight*otherCol
        fun weightedPlusAssign(otherCol: SparseIntColumn, weight: Int) {
            for(entry in otherCol.entries) {
                this.merge(entry.key, entry.value*weight) { a, b ->
                    val result = a + b
                    if(result != 0) result else null
                }
            }
        }

        fun toTrajectory(dynamics: SparseColIntMatrix, agentRows: Int, agentCols: Int): String {
            val trajectory = Array(agentCols) { Array(agentRows) { 0 } }
            for(act in this.entries) {
                val particleVal = if(act.value > 0) 1 else 2
                for(actPosition in dynamics[act.key].keys) {
                    val x = actPosition.rem(agentCols)
                    val y = actPosition.div(agentCols).rem(agentRows)
                    trajectory[x][y] = trajectory[x][y] or particleVal
                }
            }
            val particleChars = arrayOf(' ','+','-','*')
            val trajectoryString = StringBuilder()
            for(row in agentRows-1 downTo 0) {
                for(col in 0 until agentCols) {
                    trajectoryString.append(particleChars[trajectory[row][col]])
                }
                trajectoryString.appendln()
            }
            return trajectoryString.toString()
        }
    }

    class ExtendedEuclid {
        val y: Int
        val gcd: Int

        constructor(a: Int, b: Int) {
            var a = a
            var b = b
            var xp = 0
            var yp = 1
            var last_x = 1
            var last_y = 0
            var temp: Int
            while (b != 0) {
                val q = a / b
                val r = a % b
                a = b
                b = r
                temp = xp
                xp = last_x - q * xp
                last_x = temp
                temp = yp
                yp = last_y - q * yp
                last_y = temp
            }
//            x = last_x
            y = last_y
            gcd = a
        }
    }

    companion object {
        fun Identity(size: Int): SparseColIntMatrix {
            val I = SparseColIntMatrix(size,size)
            for(i in 0 until size) {
                I[i,i] = 1
            }
            return I
        }
    }
}