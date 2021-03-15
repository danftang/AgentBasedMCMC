import lib.Gnuplot
import lib.SettableLazy
import lib.gnuplot
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import kotlin.math.min

class Scratch {


    @Test
    fun stuff() {
        val set = HashSet<Double>()
        for(i in 1 until 100000) {
            var p = MinimisationMCMC.nthLogPrime(i)
            if(set.contains(p)) println("asdfasdf")
            set.add(p)

        }
        var nextFreeBasis: Int = 0
        while(!set.contains(nextFreeBasis.toDouble())) {
//                nextFreeBasis = (nextFreeBasis + 1).rem(nVariables)
            ++nextFreeBasis
            if(nextFreeBasis >= set.size) nextFreeBasis = 0
        }


    }
}