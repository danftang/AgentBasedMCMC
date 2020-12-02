package lib

open class FieldArithmetic<T>(val zero: T, val plus: (T,T)->T)

object DoubleArithmetic: FieldArithmetic<Double>(0.0, Double::plus)

object Fields {
    val doubleField = FieldArithmetic<Double>(0.0, Double::plus)
}