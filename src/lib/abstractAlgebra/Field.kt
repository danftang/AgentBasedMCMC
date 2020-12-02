package lib.abstractAlgebra

import kotlin.reflect.KClass

interface Field<T: FieldElement<T>>: org.apache.commons.math3.Field<T> {
    fun getRuntimeKClass(): KClass<T>
}