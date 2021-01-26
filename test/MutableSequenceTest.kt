import lib.collections.Removable
import lib.collections.asMutableIterable
import lib.collections.asSequenceOfRemovables
import lib.collections.mutableFlatten
import org.junit.Test

class MutableSequenceTest {

    @Test
    fun test() {
        val x = ArrayList<MutableList<Int>>()

        x.addAll(listOf(
            mutableListOf(1,2,3,4),
            mutableListOf(5),
            mutableListOf(),
            mutableListOf(6,7),
            mutableListOf()
        ))

//        val test: Sequence<Sequence<Int>>
//        test.flatten()

        val y = x.mutableFlatten()
        val it = y.iterator()
        while(it.hasNext()) {
            val element = it.next()
            println(element)
            if(element == 5) it.remove()
        }

        println(x)
    }
}