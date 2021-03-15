
fun List<Constraint<*>>.numVars(): Int {
    return this.asSequence().map { it.numVars() }.max()?:0
}

fun List<Constraint<*>>.numSlacks(): Int {
    return this.asSequence().count { it.relation != "==" }
}