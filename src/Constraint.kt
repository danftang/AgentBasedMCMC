import lib.sparseVector.SparseVector
import lib.sparseVector.asVector

interface Constraint<out COEFF> {
    val coefficients: Map<Int, COEFF>
    val relation: String
    val constant: COEFF

    fun numVars(): Int = coefficients.keys.max()?.let { it+1 }?:0



}


fun<COEFF: Comparable<COEFF>> Constraint<COEFF>.isSatisfiedBy(x: SparseVector<COEFF>): Boolean {
    val lhs = coefficients.asVector(x.operators).dotProduct(x)
    println("$lhs $relation $constant")
    return when(relation) {
        "<=" -> lhs <= constant
        ">=" -> lhs >= constant
        else -> lhs == constant
    }
}
