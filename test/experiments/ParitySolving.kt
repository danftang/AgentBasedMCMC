package experiments

import com.google.ortools.linearsolver.MPSolutionResponse
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.sat.*
import org.junit.Test
import kotlin.random.Random

class ParitySolving {

    init {
        System.loadLibrary("jniortools")
    }


    @Test
    fun testParitySolve() {
        val N = 32
        val P = 28
        val solver = CpSolver()
        val model = CpModel()
        val X = Array(N) { i ->
            model.newIntVar(0,1,"x$i")
        }

        for(j in 1..P) {
            val parity = model.newIntVar(0,N.toLong(),"parity$j")
            val randomBits = IntArray(N) { i -> Random.nextInt(0,2) }
            val constraint = model.addEquality(LinearExpr.scalProd(X + parity, randomBits + intArrayOf(-2)), 0)

            println("parity$j = ${randomBits.toList()}")

        }

        var solutions = 0
        val callback = object: CpSolverSolutionCallback() {
            override fun onSolutionCallback() {
                solutions++
                print("$solutions: ")
                for(i in 0 until N) {
                    print(value(X[i]))
                }
                println()
            }
        }
        solver.searchAllSolutions(model, callback)

    }

    @Test
    fun testParitySolveBool() {
        val N = 64
        val P = 26
        val solver = CpSolver()
        val model = CpModel()
        val X = Array(N) { i ->
            model.newBoolVar("x$i")
        }
//        model.addBoolXor(X)


        for(j in 1..P) {
            val randomBits = X.filter {_ -> Random.nextBoolean() }.toTypedArray()
            val constraint = model.addBoolXor(randomBits)
            println(randomBits.size)
        }

        var solutions = 0
        val callback = object: CpSolverSolutionCallback() {
            override fun onSolutionCallback() {
                solutions++
                print("$solutions: ")
                for(i in 0 until N) {
                    print(value(X[i]))
                }
                println()
            }
        }
        solver.searchAllSolutions(model, callback)

    }




    @Test
    fun testParitySolveMIP() {
//        val solver = MPSolver("SimpleMipProgram", MPSolver.OptimizationProblemType.CBC_MIXED_INTEGER_PROGRAMMING)
        val solver = MPSolver("SimpleMipProgram", MPSolver.OptimizationProblemType.BOP_INTEGER_PROGRAMMING)

        val N = 16
        val M = 15
        val vars = solver.makeIntVarArray(N,0.0,1.0)
        val p = solver.makeIntVarArray(M, 0.0, N.toDouble())

        for(pi in p) {
            val parity = solver.makeConstraint()
            for(v in vars) {
                if(Random.nextBoolean()) parity.setCoefficient(v,1.0)
            }
            parity.setCoefficient(pi,-2.0)
            val pval = if(Random.nextBoolean()) 1.0 else 0.0
            parity.setBounds(pval, pval)
        }

        val result = solver.solve()

        println(result)
        for(v in vars) { print("${v.solutionValue().toInt()}") }
        println()
    }

    @Test
    fun testCNF() {
        val N = 10000
        val P = 10000
        val solver = CpSolver()
        val model = CpModel()
        val X = Array(N) { i ->
            model.newBoolVar("x$i")
        }

        for(j in 1..P) {
            val randomBits = X.filter {_ -> Random.nextBoolean() }.toTypedArray()
            val constraint = model.addBoolOr(randomBits)
            println(randomBits.size)
        }

        var solutions = 0
        val callback = object: CpSolverSolutionCallback() {
            override fun onSolutionCallback() {
                solutions++
                println("$solutions: ")
//                for(i in 0 until N) {
//                    print(value(X[i]))
//                }
            }
        }
        solver.searchAllSolutions(model, callback)

    }

}