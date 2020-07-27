package lib

interface MutableMultiset<T>: MutableSet<T>, Multiset<T> {
    //override val counts: MutableMap<T, Int>

    fun retainAll(elements: Multiset<T>): Boolean

    fun compute(member: T, transform: (T,Int)->Int): Int

    fun merge(member: T, value: Int, transform: (Int,Int)->Int): Int

    override fun add(element: T) = add(element, 1)
    override fun remove(element: T) = remove(element, 1)

    fun add(m : T, n : Int) : Boolean {
        merge(m,n,Int::plus)
        return true
    }


    fun remove(element : T, n : Int) : Boolean {
        var allRemoved: Boolean = false
        merge(element, n) { oldVal, toRemove ->
            if(oldVal >= toRemove) {
                allRemoved = true
                oldVal - toRemove
            } else {
                0
            }
        }
        return allRemoved
    }


    operator fun set(m: T, n: Int) {
        compute(m) { _,_ -> n }
    }



    fun addAll(multiset: Multiset<T>): Boolean {
        multiset.counts.forEach { add(it.key, it.value) }
        return true
    }

    override fun addAll(elements: Collection<T>): Boolean {
        var success = true
        elements.forEach {
            success = success && add(it)

        }
        return success
    }

    override fun removeAll(elements: Collection<T>): Boolean {
        var success = true
        elements.forEach {
            success = success && remove(it)

        }
        return success
    }

    override fun retainAll(elements: Collection<T>): Boolean {
        return retainAll(elements.toMutableMultiset())
    }

}

fun <T>Collection<T>.toMutableMultiset(): MutableMultiset<T> {
    return HashMultiset(this)
}

fun <T>Iterable<T>.toMutableMultiset(): MutableMultiset<T> {
    return HashMultiset(this)
}

fun <T>Sequence<T>.toMutableMultiset(): MutableMultiset<T> {
    return HashMultiset(this)
}

fun <T>Array<out T>.toMutableMultiset(): MutableMultiset<T> {
    return HashMultiset(this)
}
