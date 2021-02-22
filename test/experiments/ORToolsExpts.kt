package experiments

import Constraint
import org.junit.Test

class ORToolsExpts {
    @Test
    fun solveTest() {
        val result = ORTools.GlopSolve(
            listOf(
                Constraint(mutableMapOf(0 to 1),"<=", 1),
                Constraint(mutableMapOf(0 to 1, 1 to 1),"<=", 2)
            ),
            mapOf(1 to -1)
        )
        println(result.asList())

        val intResult = ORTools.IntegerSolve(
            listOf(
                Constraint(mutableMapOf(0 to 1),"<=", 1),
                Constraint(mutableMapOf(0 to 1, 1 to 1),"<=", 2)
            ),
            mapOf(1 to -1)
        )
        println(intResult.asList())

    }
}