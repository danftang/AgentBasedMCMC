import lib.collections.Multiset
import org.junit.Test

class TrajectoryTest {
    @Test
    fun testPairs() {
        val t = Trajectory<CatAndMouseABM.CatMouseAgent,CatAndMouseABM.Acts>()

//        t[0,CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT),CatAndMouseABM.Acts.STAYPUT] = 1

        val pair = Pair(CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT),CatAndMouseABM.Acts.STAYPUT)
        t.add(Multiset())
        t[0][pair] = 1

        println(t)
        println(t[0,CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT),CatAndMouseABM.Acts.STAYPUT])
        println(t[0].nonZeroEntries[Pair(CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT),CatAndMouseABM.Acts.STAYPUT)])
        println(t[0][pair])

        println(pair.second == CatAndMouseABM.Acts.STAYPUT)
        println(pair.first == CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT))
        println(Pair(CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT),CatAndMouseABM.Acts.STAYPUT) == pair)
    }
}