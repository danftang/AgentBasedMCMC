
class Predator : Agent {
    constructor(pos: Int) : super(pos)


    override fun copyAt(pos: Int) = Predator(pos)


    fun hamiltonian(h: Hamiltonian<Agent>, params: Params) {
        if(params.predDiffuse > 0.0) diffuse(h, params.GRIDSIZE, params.predDiffuse)
        if(params.predDie > 0.0) die(h, params.predDie)
//        if(params.predCaptureOnly > 0.0) capture(h, params)
//        if(params.predCaptureAndReproduce > 0.0) captureAndReproduce(h, params)
    }

    fun diffuse(h: Hamiltonian<Agent>, size: Int, rate: Double) {
        h += action(rate/5.0, copyAt(right(size)))
        h += action(rate/5.0, copyAt(left(size)))
        h += action(rate/5.0, copyAt(up(size)))
        h += action(rate/5.0, copyAt(down(size)))
        h += action(rate/5.0, this)
    }

//    fun capture(h: Hamiltonian<Agent>, params: Params) {
//        h += interaction(params.predCaptureOnly, Prey(pos), this)
//        h += interaction(params.predCaptureOnly, Prey(right(params.GRIDSIZE)), this)
//        h += interaction(params.predCaptureOnly, Prey(left(params.GRIDSIZE)), this)
//        h += interaction(params.predCaptureOnly, Prey(up(params.GRIDSIZE)), this)
//        h += interaction(params.predCaptureOnly, Prey(down(params.GRIDSIZE)), this)
//    }
//
//    fun captureAndReproduce(h: Hamiltonian<Agent>, params: Params) {
//        h += interaction(params.predCaptureAndReproduce, Prey(right(params.GRIDSIZE)), this, Predator(right(params.GRIDSIZE)))
//        h += interaction(params.predCaptureAndReproduce, Prey(left(params.GRIDSIZE)), this, Predator(left(params.GRIDSIZE)))
//        h += interaction(params.predCaptureAndReproduce, Prey(up(params.GRIDSIZE)), this, Predator(up(params.GRIDSIZE)))
//        h += interaction(params.predCaptureAndReproduce, Prey(down(params.GRIDSIZE)), this, Predator(down(params.GRIDSIZE)))
//    }

    override fun toString() = "f($pos)"

    override fun hashCode() = pos*2 + 1

    override fun equals(other: Any?) = ((other is Predator) && (pos == other.pos))
}