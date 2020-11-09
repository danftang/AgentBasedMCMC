import org.junit.Test

class Scratch {


    data class MatrixEntry(var value: Int)

//    class ColEntry(val col: Int, row: Int, value: Int): RowEntry(row, value) {
//        override fun hashCode(): Int {
//            return col
//        }
//    }

    @Test
    fun stuff() {
        val row = HashMap<Int,MatrixEntry>()
        val col = HashMap<Int,MatrixEntry>()


        val newEntry = MatrixEntry(1234)
        col[2] = newEntry
        row[1] = newEntry

        col[2]?.value = 2345
//        row[1]?.value = 5432

        println(col[2]?.value)
        println(row[1]?.value)

    }
}