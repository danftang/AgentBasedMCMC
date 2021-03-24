
fun<T: Number> List<MutableConstraint<T>>.numVars(): Int {
    return this.asSequence().map { it.numVars() }.max()?:0
}

fun<T: Number> List<MutableConstraint<T>>.numSlacks(): Int {
    return this.count { it.relation != "==" }
}