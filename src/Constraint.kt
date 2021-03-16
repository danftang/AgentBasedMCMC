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


// returns the absolute value of the difference between lhs and rhs
//fun<COEFF: Number> Constraint<COEFF>.slackness(values: SparseVector<COEFF>): Double {
//    if(relation == "==") return 0.0
//    val lhs = coefficients.asVector(values.operators).dotProduct(values).toDouble()
//    return with(values.operators) {
//        if (relation == "<=") constant.toDouble() - lhs else lhs - constant.toDouble()
//    }
//}

// returns the absolute value of the difference between lhs and rhs
fun<COEFF: Number> Constraint<COEFF>.slackness(values: Map<Int,Double>): Double {
    if(relation == "==") return 0.0
    val lhs = values.entries.sumByDouble { coefficients[it.key]?.toDouble()?:0.0 * it.value }
    return if (relation == "<=") constant.toDouble() - lhs else lhs - constant.toDouble()
}
