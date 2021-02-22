
fun<T: Number> List<Constraint<T>>.numVars(): Int {
    return this.asSequence().map { it.numVars() }.max()?:0
}