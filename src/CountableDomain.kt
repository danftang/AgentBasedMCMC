interface CountableDomain<T> {
    val size: Int

    operator fun get(index: Int): T
}

inline fun<reified T : Enum<T>> enumCountableDomain() = object: CountableDomain<T> {
    override val size: Int
        get() = enumValues<T>().size
//    override fun toIndex(agent: T) = agent.ordinal
    override fun get(index: Int) = enumValues<T>()[index]
}
