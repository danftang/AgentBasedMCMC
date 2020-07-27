package cnf

class Literal(val id: Int) {

    operator fun not() = Literal(-id)

    operator fun plus(other: Literal) = OrClause(this, other)

    companion object {
        var idCounter: Int = 0

        val nextId: Int
                get() = ++idCounter
    }

    override fun toString(): String {
        return id.toString()
    }
}