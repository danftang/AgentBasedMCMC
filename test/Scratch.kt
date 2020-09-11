import org.junit.Test

class Scratch {

    data class MatrixEntry(var value: Int)

    @Test
    fun stuff() {
        val map1 = HashMap<Int,MatrixEntry>()
        val map2 = HashMap<Int,MatrixEntry>()
        val myEntry = MatrixEntry(1234)
        map1[1] = myEntry
        map2[2] = myEntry


        println("${map1[1]} ${map2[2]}")
        map1[1]?.value = 2345
        println("${map1[1]} ${map2[2]}")

    }
}