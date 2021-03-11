import lib.Gnuplot
import lib.SettableLazy
import lib.gnuplot
import org.junit.Test
import kotlin.math.min

class Scratch {


    @Test
    fun stuff() {
        var my: Double by SettableLazy { 0.1234 }
//        my = 2.345
        println(my)
    }
}