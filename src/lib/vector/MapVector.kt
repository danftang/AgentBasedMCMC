package lib.vector

import lib.abstractAlgebra.FieldOperators

class MapVector<T: Any>(
    val elementOperators: FieldOperators<T>,
    override val nonZeroEntries: Map<Int,T>)
    : SparseVector<T>, FieldOperators<T> by elementOperators
{
    override fun new() = MutableMapVector(elementOperators)
}
