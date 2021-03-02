package lib.sparseVector

import lib.abstractAlgebra.FieldOperators

class MutableMapVector<T>(
    val fieldOperators: FieldOperators<T>,
    override val nonZeroEntries: MutableMap<Int,T> = HashMap())
    : MutableSparseVector<T>, FieldOperators<T> by fieldOperators
{

    constructor(copyFrom: SparseVector<T>): this(copyFrom.operators, HashMap(copyFrom.nonZeroEntries))

    override fun new() = MutableMapVector(fieldOperators)
    override fun toString() = nonZeroEntries.toString()
    override fun equals(other: Any?) = isEqualTo(other)
}

inline fun<T: Any> MutableMap<Int,T>.asMutableVector(operators: FieldOperators<T>): MutableMapVector<T> =
    MutableMapVector(operators, this)

//inline fun<T: FieldElement<T>> MutableMap<Int,T>.asMutableMapVector(apacheField: Field<T>): MutableMapVector<T> =
//    MutableMapVector(apacheField.asFieldOperators(), this)
