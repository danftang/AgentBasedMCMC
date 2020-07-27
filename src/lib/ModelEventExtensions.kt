package lib

import Event

fun <T> Multiset<Event<T>>.consequences(): Multiset<T> {
    return this.asSequence().flatMap { it.consequences.asSequence() }.toMultiset()
}

fun <T> Multiset<Event<T>>.primaryRequirements(): Multiset<T> {
    return this.asSequence().map { it.primaryAgent }.toMultiset()
}

fun <T> Set<Event<T>>.consequences(): Multiset<T> {
    return this.asSequence().flatMap { it.consequences.asSequence() }.toMultiset()
}

fun <T> Set<Event<T>>.primaryRequirements(): Multiset<T> {
    return this.asSequence().map { it.primaryAgent }.toMultiset()
}
