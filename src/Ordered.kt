interface Ordered<T> {
    val ordinal: Int
    val domain: CountableDomain<T>
}