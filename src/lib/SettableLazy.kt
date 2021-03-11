package lib

import kotlin.properties.ReadWriteProperty
import kotlin.reflect.KProperty

class SettableLazy<OBJ,T>(val initialiser: ()->T): ReadWriteProperty<OBJ,T> {
    var backing: T? = null

    override fun getValue(thisRef: OBJ, property: KProperty<*>): T {
        if(backing == null) { backing = initialiser() }
        return backing!! // NOT thread safe but hey ho
    }

    override fun setValue(thisRef: OBJ, property: KProperty<*>, value: T) {
        backing = value
    }
}