package experiments

import com.google.ortools.sat.*
import org.junit.Test
import kotlin.random.Random

class ProbabilityStratification {

    init { System.loadLibrary("jniortools") }

    @Test
    fun testABMtoSat() {
        val N = 16
        val T = 3

        val solver = CpSolver()
        val model = CpModel()
        val logprobs = HashMap<IntVar,Int>()

//        var T0 = Array(N) { Array(N) { context.newVar() } }
//        var T1: Array<Array<Literal>>

        val occupation = ArrayList<Array<Array<IntVar>>>()

        occupation.add(Array(N) { Array(N) { i -> model.newIntVar(0,4,"0$i") } })

        for(t in 1..T) {
            occupation.add(Array(N) { Array(N) { i -> model.newIntVar(0,4,"$t$i")} })
//            T1 = Array(N) { Array(N) { context.newDependentVar() } }

            val actsByStartState = Array(N) { Array(N) { ArrayList<IntVar>() } }
            val actsByEndState = Array(N) { Array(N) { ArrayList<IntVar>() } }

            for (row in 0 until N) {
                for (col in 0 until N) {
                    val left = model.newBoolVar("l$t$row$col")
                    val right = model.newBoolVar("r$t$row$col")
                    val up = model.newBoolVar("u$t$row$col")
                    val down = model.newBoolVar("d$t$row$col")
                    actsByStartState[row][col].addAll(listOf(left, right, up, down))
                    actsByEndState[row][(col + N - 1).rem(N)].add(left)
                    actsByEndState[row][(col + 1).rem(N)].add(right)
                    actsByEndState[(row + 1).rem(N)][col].add(up)
                    actsByEndState[(row + N - 1).rem(N)][col].add(down)
                    logprobs[left] = Random.nextInt(10)
                    logprobs[right] = Random.nextInt(10)
                    logprobs[up] = Random.nextInt(10)
                    logprobs[down] = Random.nextInt(10)
                }
            }

            for (row in 0 until N) {
                for (col in 0 until N) {
                    model.addEquality(
                        LinearExpr.sum(actsByStartState[row][col].toTypedArray()),
                        LinearExpr.term(occupation[t-1][row][col],1)
                    )
                    model.addEquality(
                        LinearExpr.sum(actsByEndState[row][col].toTypedArray()),
                        LinearExpr.term(occupation[t][row][col],1)
                    )
                }
            }

        }

        // observations
        for (row in 0 until N) {
            for (col in 0 until N) {
                model.addEquality(
                    LinearExpr.term(occupation.last()[row][col],1),
                    if(Random.nextDouble() < 0.1) 1L else 0L
                )
            }
        }

        // ProbStrat

        model.addLinearConstraint(LinearExpr.scalProd(logprobs.keys.toTypedArray(), logprobs.values.toIntArray()), 200, 200)

        var solutions = 0
        val callback = object: CpSolverSolutionCallback() {
            override fun onSolutionCallback() {
                solutions++
                print("$solutions: ")
//                for(i in 0 until N) {
//                    print(value(X[i]))
//                }
                println()
            }
        }

        println(model.modelStats())
        solver.searchAllSolutions(model, callback)

    }
}