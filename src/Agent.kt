import lib.collections.Multiset
import kotlin.random.Random

interface  Agent<AGENT>: Ordered<AGENT> {
    fun timestep(others: Multiset<AGENT>): Array<Double>
}