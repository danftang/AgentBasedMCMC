package lib.abstractAlgebra

inline class DoubleField(val x: Double): FieldElement<DoubleField>, Comparable<DoubleField> {

    override fun plus(other: DoubleField)  = DoubleField(x + other.x)
    override fun minus(other: DoubleField) = DoubleField(x - other.x)
    override fun unaryMinus()              = DoubleField(-x)
    override fun times(n: Int)             = DoubleField(x * n)
    override fun times(other: DoubleField) = DoubleField(x * other.x)
    override fun div(other: DoubleField)   = DoubleField(x / other.x)
    override fun getField()                = DoubleField
    override fun compareTo(other: DoubleField) = x.compareTo(other.x)

    fun asDouble(): Double = x

    override fun toString() = x.toString()

    companion object: Field<DoubleField> {
        override fun getZero()          = DoubleField(0.0)
        override fun getOne()           = DoubleField(1.0)
        override fun getRuntimeClass()  = DoubleField::class.java
        override fun getRuntimeKClass() = DoubleField::class

        inline fun plus(a: Double, b: Double): Double = a + b // test
    }

}
