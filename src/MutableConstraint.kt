import lib.sparseVector.SparseVector
import lib.sparseVector.asVector
import java.lang.IllegalArgumentException
import java.lang.StringBuilder

class MutableConstraint<COEFF>: Constraint<COEFF> {
    override val coefficients: MutableMap<Int, COEFF>
    override val relation: String
    override val constant: COEFF

    constructor(coefficients: MutableMap<Int, COEFF>, relation: String, constant: COEFF) {
        this.coefficients = coefficients
        this.relation = when(relation) {
            "==", "<=", ">=" -> relation
            else -> throw(IllegalArgumentException("Relation should be one of '==','<=' or '>='"))
        }
        this.constant = constant
    }

    // returns the absolute velue of the difference between lhs and rhs
    fun slackness(values: SparseVector<COEFF>): COEFF {
        if(relation == "==") return values.operators.zero
        val lhs = coefficients.asVector(values.operators).dotProduct(values)
        return with(values.operators) {
            if (relation == "<=") constant - lhs else lhs - constant
        }
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

