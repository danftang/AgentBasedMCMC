package lib.abstractAlgebra

inline class IntField(val x: Int): FieldElement<IntField> {

    override fun plus(other: IntField)  = IntField(x + other.x)
    override fun minus(other: IntField) = IntField(x - other.x)
    override fun unaryMinus()           = IntField(-x)
    override fun times(n: Int)          = IntField(x * n)
    override fun times(other: IntField) = IntField(x * other.x)
    override fun div(other: IntField)   = IntField(x / other.x)
    override fun getField()             = IntField

    override fun compareTo(other: IntField) = x.compareTo(other.x)

    fun asInt(): Int = x

    companion object: Field<IntField> {
        override fun getZero()          = IntField(0)
        override fun getOne()           = IntField(1)
        override fun getRuntimeClass()  = IntField::class.java
        override fun getRuntimeKClass() = IntField::class
    }

}