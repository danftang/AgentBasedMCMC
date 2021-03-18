import lib.Gnuplot
import lib.SettableLazy
import lib.gnuplot
import org.junit.Test
import java.time.Instant
import kotlin.math.min
import kotlin.math.sin
import kotlin.math.sqrt

class Scratch {


    @Test
    fun stuff() {
        val inst = Instant.now()
        println(inst.toEpochMilli())
        val d = Array(10000000) { sin(sqrt(it.toDouble())) }
        println(Instant.now().toEpochMilli())
    }
}