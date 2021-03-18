package experiments

import MutableConstraint
import org.junit.Test

class ORToolsExpts {
    @Test
    fun solveTest() {
        val result = ORTools.GlopSolve(
            listOf(
                MutableConstraint(mutableMapOf(0 to 1),"<=", 1),
                MutableConstraint(mutableMapOf(0 to 1, 1 to 1),"<=", 2)
            ),
            mapOf(1 to -1)
        )
        println(result.asList())

        val intResult = ORTools.IntegerSolve(
            listOf(
                MutableConstraint(mutableMapOf(0 to 1),"<=", 1),
                MutableConstraint(mutableMapOf(0 to 1, 1 to 1),"<=", 2)
            ),
            mapOf(1 to -1)
        )
        println(intResult.asList())

    }
}