package lib.collections

class Multiset<T>(val entries: MutableMap<T,Int> = HashMap()) { //: MutableIterable<MutableMap.MutableEntry<T,Int>> {

    val members: MutableSet<T>
        get() = entries.keys

    val occupationNumbers: MutableCollection<Int>
        get() = entries.values


    operator fun set(agent: T, occupation: Int) {
        if(occupation != 0) {
            entries[agent] = occupation
        } else {
            entries.remove(agent)
        }
    }

    operator fun get(agent: T): Int {
        return entries.getOrDefault(agent, 0)
    }

    operator fun plusAssign(other: Multiset<T>) {
        other.entries.forEach {
            this[it.key] += it.value
        }
    }

//    override fun iterator(): MutableIterator<MutableMap.MutableEntry<T, Int>> {
//        return nonZeroEntries.iterator()
//    }

    override fun toString(): String {
        return entries.toString()
    }

}

fun<T> multisetOf(vararg nonZeroEntries: Pair<T,Int>): Multiset<T> = Multiset(hashMapOf(*nonZeroEntries))