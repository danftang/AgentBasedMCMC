
fun<T: Number> List<MutableConstraint<T>>.numVars(): Int {
    return this.asSequence().map { it.numVars() }.max()?:0
}