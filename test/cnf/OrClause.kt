package cnf

import java.lang.StringBuilder
import kotlin.math.absoluteValue

class OrClause: ArrayList<Literal> {

    constructor(vararg literals: Literal) {
        addAll(literals)
    }

    constructor(literals: Collection<Literal>) {
        addAll(literals)
    }

    constructor(literals: Sequence<Literal>) {
        addAll(literals)
    }

    fun nMax(): Int {
        var max = 0
        forEach {
            val v = it.id.absoluteValue
            if(v > max) max = v
        }
        return max
    }

    operator fun plus(v: Literal): OrClause {
        val newClause = OrClause()
        newClause.addAll(this)
        newClause.add(v)
        return newClause
    }

    operator fun times(other: OrClause): CnfFormula {
        val newFormula = CnfFormula()
        newFormula.add(this)
        newFormula.add(other)
        return(newFormula)
    }

    override fun toString(): String {
        val s = StringBuilder()
        this.forEach { lit ->
            s.append(lit.toString())
            s.append(' ')
        }
        s.append('0')
        return s.toString()
    }
}