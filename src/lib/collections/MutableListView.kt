package lib.collections

class MutableListView<STORE,VIEW>(val underlyingList: MutableList<STORE>, val readTransform: (STORE)->VIEW): AbstractMutableList<VIEW>() {
    override val size: Int
        get() = underlyingList.size

    override fun add(index: Int, element: VIEW) {
        TODO("Not yet implemented")
    }

    override fun get(index: Int): VIEW {
        return readTransform(underlyingList[index])
    }

    override fun removeAt(index: Int): VIEW {
        return readTransform(underlyingList.removeAt(index))
    }

    override fun set(index: Int, element: VIEW): VIEW {
        TODO("Not yet implemented")
    }

}