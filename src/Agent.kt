import lib.collections.Multiset
import kotlin.random.Random

interface  Agent<AGENT> {
    fun timestep(others: Multiset<AGENT>): Array<Double>
}