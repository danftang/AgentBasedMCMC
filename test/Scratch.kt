import lib.Gnuplot
import lib.gnuplot
import org.junit.Test

class Scratch {


    @Test
    fun stuff() {
//        gnuplot {
//            invoke("plot x*x with lines")
//        }

        val gnuplot = Gnuplot()
        with(gnuplot) {
            invoke("plot x*x with lines")
            close()
        }

    }
}