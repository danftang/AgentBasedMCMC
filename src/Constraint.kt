import java.lang.IllegalArgumentException
import java.lang.StringBuilder

class Constraint<COEFF> {
    val coefficients: MutableMap<Int, COEFF>
    val relation: String
    val constant: COEFF

    constructor(coefficients: MutableMap<Int, COEFF>, relation: String, constant: COEFF) {
        this.coefficients = coefficients
        this.relation = when(relation) {
            "==", "<=", ">=" -> relation
            else -> throw(IllegalArgumentException("Relation should be one of '==','<=' or '>='"))
        }
        this.constant = constant
    }

    override fun toString(): String {
        val out = StringBuilder()
        for(entry in coefficients) {
            out.append(" ${entry.value}x_${entry.key} +")
        }
        if(out.lastIndex > 0) out.deleteCharAt(out.lastIndex)
        out.append(" $relation ")
        out.append(constant.toString())
        return out.toString()
    }
}