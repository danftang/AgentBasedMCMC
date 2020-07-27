package cnf

import java.lang.StringBuilder

class CnfContext {
    private var lastVarId = 0
    val independentVars = ArrayList<Literal>()

    val nVars: Int
        get() = lastVarId

    fun newVar(): Literal {
        val v = Literal(++lastVarId)
        independentVars.add(v)
        return v
    }

    fun newDependentVar() = Literal(++lastVarId)


    // Create a new dependent variable such that x = a or b or c or...
    // where literals = (a,b,c...)
    // returns x and a formula that asserts the equality.
    //
    // In CNF, this translates to
    // (!x + a + b + c + ...)*(x + !a)*(x + !b)*(x + !c)...
    fun newDependentOr(vararg literals: Literal): Pair<Literal,CnfFormula> {
        val x = newDependentVar()
        val f = CnfFormula()
        val notxClause = OrClause(literals.asSequence() + !x)
        f.add(notxClause)
        for(v in literals) {
            f.add(x + !v)
        }
        return Pair(x, f)
    }

    override fun toString(): String {
        val s = StringBuilder()
        s.append("c ind ")
        for(v in independentVars) {
            s.append(v)
            s.append(' ')
        }
        s.append('0')
        return(s.toString())
    }
}