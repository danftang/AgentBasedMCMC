interface CountableDomain<T> {
    val size: Int

    fun toIndex(agent: T): Int
    fun toObject(index: Int): T

}

inline fun<reified T : Enum<T>> countableDomainOf() = object: CountableDomain<T> {
    override val size: Int
        get() = enumValues<T>().size
    override fun toIndex(agent: T) = agent.ordinal
    override fun toObject(index: Int) = enumValues<T>()[index]
}
