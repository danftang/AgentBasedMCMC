interface Tableaux<T> {
    val rows: List<Map<Int,T>>
    val columns: List<Map<Int,T>>

    val nCols: Int
    val nRows: Int

    operator fun get(i: Int, j: Int): T
    operator fun set(i: Int, j: Int)
}