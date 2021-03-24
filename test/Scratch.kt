import lib.Gnuplot
import lib.SettableLazy
import lib.gnuplot
import org.apache.commons.math3.linear.OpenMapRealVector
import org.apache.commons.math3.linear.RealVectorChangingVisitor
import org.apache.commons.math3.util.OpenIntToDoubleHashMap
import org.junit.Test
import java.time.Instant
import kotlin.math.min
import kotlin.math.sin
import kotlin.math.sqrt

class Scratch {


    @Test
    fun stuff() {
        val m = OpenIntToDoubleHashMap(0.0)
        m.put(1, 1234.0)
        m.put(10,456.0)
        m.put(5,4567.0)
        m.put(2,0.0)
        println(m)

        var it = m.iterator()
        while(it.hasNext()) {
            it.advance()
            println("${it.key()} ${it.value()}")
        }

//        it = m.iterator()
//        while(it.hasNext()) {
//            it.advance()
//            m.put(it.key(), it.value()*0.5)
//        }

//        val iter = m.iterator()
//        val keys = IntArray(m.size()) {
//            iter.advance()
//            iter.key()
//        }
//        keys.forEach { m.remove(it) }
//
//        println()
//        it = m.iterator()
//        while(it.hasNext()) {
//            it.advance()
//            println("${it.key()} ${it.value()}")
//        }
//

//        val v = OpenMapRealVector(15)
//        v.setEntry(1, 1234.0)
//        v.setEntry(10, 456.0)
//
//        val visitor = object : RealVectorChangingVisitor {
//            override fun start(p0: Int, p1: Int, p2: Int) {
//            }
//
//            override fun visit(p0: Int, p1: Double): Double {
//                return p1 * 0.5
//            }
//
//            override fun end(): Double {
//                return 1.0
//            }
//
//        }
//        for(i in 0 until 15) print("${v.getEntry(i)} ")
//        println()
//
//        v.walkInOptimizedOrder(visitor)
//
//        for(i in 0 until 15) print("${v.getEntry(i)} ")
//        println()
//
//        v.mapToSelf { it * 0.5 }
//
//        for(i in 0 until 15) print("${v.getEntry(i)} ")
//        println()
    }
}