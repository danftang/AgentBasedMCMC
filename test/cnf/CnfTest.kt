package cnf

import org.junit.Test

class CnfTest {

   @Test
    fun testCnf() {
       val context = CnfContext()
       val x = context.newVar()
       val y = context.newVar()
       val z = context.newVar()
       val a = context.newVar()

       val f = CnfFormula()
       f.addDemux(context, a,x,y,z)
       val satProblem = CnfSatProblem(context, f)

       println(satProblem.toString())
    }

    @Test
    fun testABMtoCnf() {
        val N = 8
        val T = 8
        val context = CnfContext()
        val f = CnfFormula()

//        var T0 = Array(N) { Array(N) { context.newVar() } }
//        var T1: Array<Array<Literal>>

        val occupation = ArrayList<Array<Array<Literal>>>()

        occupation.add(Array(N) { Array(N) { context.newDependentVar() } })

        for(t in 1..T) {
            occupation.add(Array(N) { Array(N) { context.newDependentVar() } })
//            T1 = Array(N) { Array(N) { context.newDependentVar() } }

            val actsByStartState = Array(N) { Array(N) { ArrayList<Literal>() } }
            val actsByEndState = Array(N) { Array(N) { ArrayList<Literal>() } }

            for (row in 0 until N) {
                for (col in 0 until N) {
                    val left = context.newVar()
                    val right = context.newVar()
                    val up = context.newVar()
                    val down = context.newVar()
                    actsByStartState[row][col].addAll(listOf(left, right, up, down))
                    actsByEndState[row][(col + N - 1).rem(N)].add(left)
                    actsByEndState[row][(col + 1).rem(N)].add(right)
                    actsByEndState[(row + 1).rem(N)][col].add(up)
                    actsByEndState[(row + N - 1).rem(N)][col].add(down)
                }
            }

            for (row in 0 until N) {
                for (col in 0 until N) {
                    f.addDemux(context, occupation[t-1][row][col], *actsByStartState[row][col].toTypedArray())
                    f.addOr(occupation[t][row][col], *actsByEndState[row][col].toTypedArray())
                }
            }

        }

        // observations
        for (row in 0 until N) {
            for (col in 0 until N) {
                if(row == N/2 && col ==N/2) {
                    f.add(OrClause(occupation.last()[row][col]))
                    f.add(OrClause(occupation[0][row][col]))
                } else {
                    f.add(OrClause(!occupation.last()[row][col]))
                    f.add(OrClause(!occupation[0][row][col]))
                }
            }
        }
        val problem = CnfSatProblem(context, f)
        println(problem)
        problem.save("abm.cnf")
    }
}