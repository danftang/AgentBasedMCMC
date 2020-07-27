package cnf

import java.io.File
import java.io.FileWriter
import java.lang.StringBuilder

class CnfSatProblem(val context: CnfContext, val formula: CnfFormula) {

    fun save(filename: String) {
        val file = FileWriter(filename)
        file.write(toString())
        file.close()
    }

    override fun toString(): String {
        val s = StringBuilder()
        s.appendln(context.toString())
        s.appendln("p cnf ${context.nVars} ${formula.size}")
        s.append(formula.toString())
        return s.toString()
    }
}