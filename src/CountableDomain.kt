interface CountableDomain<T>: Iterable<T> {
    val size: Int

    operator fun get(index: Int): T

    override fun iterator(): Iterator<T> {
        return object: Iterator<T> {
            var index: Int = 0
            override fun hasNext() = index < size
            override fun next() = get(index++)
        }
    }
}

inline fun<reified T : Enum<T>> enumCountableDomain() = object: CountableDomain<T> {
    override val size: Int
        get() = enumValues<T>().size
//    override fun toIndex(agent: T) = agent.ordinal
    override fun get(index: Int) = enumValues<T>()[index]
}
