package lib.vector

import lib.abstractAlgebra.FieldOperators

class MutableMapVector<T: Any>(
    val elementOperators: FieldOperators<T>,
    override val nonZeroEntries: MutableMap<Int,T> = HashMap())
    : MutableSparseVector<T>, FieldOperators<T> by elementOperators
{
    override fun new() = MutableMapVector(elementOperators)
}
