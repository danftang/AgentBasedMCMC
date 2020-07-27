package cnf

import java.lang.StringBuilder
import kotlin.math.absoluteValue

class CnfFormula: ArrayList<OrClause> {

    constructor(vararg clauses: OrClause) {
        addAll(clauses)
    }

    constructor(clauses: Collection<OrClause>) {
        addAll(clauses)
    }

    // Add a demultiplexer between input and outputs
    // i.e. if i = input and (a,b,c..) = outputs
    // i = (a + b + c)
    // no more than one of a,b,c... is true
    fun addDemux(context: CnfContext, input: Literal, vararg outputs: Literal) {
        if(outputs.size == 2) {
            addDemux(input, outputs[0], outputs[1])
            return
        }
        val reducedOutputs = Array((outputs.size+1)/2) { index ->
            if(2*index + 1 >= outputs.size) {
                outputs[2*index]
            } else {
                val intermediateDeMux = context.newDependentVar()
                addDemux(intermediateDeMux, outputs[2 * index], outputs[2 * index + 1])
                intermediateDeMux
            }
        }
        addDemux(context, input, *reducedOutputs)
    }


    fun addDemux(input: Literal, output1: Literal, output2: Literal) {
        addOr(input, output1, output2)
        add(!output1 + !output2)
    }


    // Assert that x = a or b or c or...
    // where result = x and literals = (a,b,c...)
    //
    // In CNF, this translates to
    // (!x + a + b + c + ...)*(x + !a)*(x + !b)*(x + !c)...
    fun addOr(result: Literal, vararg inputs: Literal) {
        val notxClause = OrClause(inputs.asSequence() + !result)
        add(notxClause)
        for(v in inputs) {
            add(result + !v)
        }
    }


    operator fun times(clause: OrClause): CnfFormula {
        val newFormula = CnfFormula()
        newFormula.addAll(this)
        newFormula.add(clause)
        return newFormula
    }




    override fun toString(): String {
        val s = StringBuilder()
        for(clause in this) {
            s.appendln(clause.toString())
        }
        return s.toString()
    }
}