package lib

import java.io.*

fun File.objectOutputStream(): ObjectOutputStream {
    return ObjectOutputStream(this.outputStream())
}

fun File.objectInputStream(): ObjectInputStream {
    return ObjectInputStream(this.inputStream())
}

fun File.writeObject(obj: Any) {
    this.objectOutputStream().use { it.writeObject(obj) }
}

fun <T>File.readObject(): T {
    return this.objectInputStream().use { it.readObject() as T }
}