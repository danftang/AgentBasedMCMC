package lib.vector

import lib.abstractAlgebra.DoubleOperators
import lib.abstractAlgebra.FieldOperators
import lib.abstractAlgebra.IntOperators
import lib.abstractAlgebra.asFieldOperators
import org.apache.commons.math3.Field
import org.apache.commons.math3.FieldElement

class MutableMapVector<T>(
    val fieldOperators: FieldOperators<T>,
    override val nonZeroEntries: MutableMap<Int,T> = HashMap())
    : MutableSparseVector<T>, FieldOperators<T> by fieldOperators
{
    override fun new() = MutableMapVector(fieldOperators)
    override fun toString() = nonZeroEntries.toString()
    override fun equals(other: Any?) = isEqualTo(other)
}

inline fun<T: Any> MutableMap<Int,T>.asMutableVector(operators: FieldOperators<T>): MutableMapVector<T> =
    MutableMapVector(operators, this)

//inline fun<T: FieldElement<T>> MutableMap<Int,T>.asMutableMapVector(apacheField: Field<T>): MutableMapVector<T> =
//    MutableMapVector(apacheField.asFieldOperators(), this)
