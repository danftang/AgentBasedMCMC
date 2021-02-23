package lib.collections

class Multiset<T>(val nonZeroEntries: MutableMap<T,Int> = HashMap()): MutableIterable<MutableMap.MutableEntry<T,Int>> {

    val members: MutableSet<T>
        get() = nonZeroEntries.keys

    val occupationNumbers: MutableCollection<Int>
        get() = nonZeroEntries.values

    operator fun set(agent: T, occupation: Int) {
        if(occupation != 0) {
            nonZeroEntries[agent] = occupation
        } else {
            nonZeroEntries.remove(agent)
        }
    }

    operator fun get(agent: T): Int {
        return nonZeroEntries.getOrDefault(agent, 0)
    }

    override fun iterator(): MutableIterator<MutableMap.MutableEntry<T, Int>> {
        return nonZeroEntries.iterator()
    }


}