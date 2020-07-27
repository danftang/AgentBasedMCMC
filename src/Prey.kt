
class Prey : Agent {
    constructor(pos: Int) : super(pos)

    override fun copyAt(pos: Int) = Prey(pos)

    fun hamiltonian(h: Hamiltonian<Agent>, params: Params) {
        reproduce(h, params)
        diffuse(h, params)
        die(h, params.preyDie)
        beEaten(h, params)
        beEatenAndReplaced(h, params)
    }

    fun diffuse(h: Hamiltonian<Agent>, params: Params) {
        val rate = params.preyDiffuse
        val size = params.GRIDSIZE
        h += action(rate/5.0, setOf(Predator(right(size))), Prey(right(size)))
        h += action(rate/5.0, setOf(Predator(left(size))), Prey(left(size)))
        h += action(rate/5.0, setOf(Predator(up(size))), Prey(up(size)))
        h += action(rate/5.0, setOf(Predator(down(size))), Prey(down(size)))
        h += action(rate/5.0, setOf(Predator(pos)), this)
    }

    fun reproduce(h: Hamiltonian<Agent>, params: Params) {
        if(params.preyReproduce > 0.0) {
            h += action(params.preyReproduce / 4.0, Prey(left(params.GRIDSIZE)), this)
            h += action(params.preyReproduce / 4.0, Prey(right(params.GRIDSIZE)), this)
            h += action(params.preyReproduce / 4.0, Prey(up(params.GRIDSIZE)), this)
            h += action(params.preyReproduce / 4.0, Prey(down(params.GRIDSIZE)), this)
        }
    }

    fun beEaten(h: Hamiltonian<Agent>, params: Params) {
        if(params.predCaptureOnly > 0.0) {
            h += interaction(params.predCaptureOnly / 5.0, Predator(right(params.GRIDSIZE)))
            h += interaction(params.predCaptureOnly / 5.0, Predator(left(params.GRIDSIZE)))
            h += interaction(params.predCaptureOnly / 5.0, Predator(up(params.GRIDSIZE)))
            h += interaction(params.predCaptureOnly / 5.0, Predator(down(params.GRIDSIZE)))
            h += interaction(params.predCaptureOnly / 5.0, Predator(pos))
        }
    }

    fun beEatenAndReplaced(h: Hamiltonian<Agent>, params: Params) {
        if(params.predCaptureAndReproduce > 0.0) {
            h += interaction(params.predCaptureAndReproduce / 5.0, Predator(right(params.GRIDSIZE)), Predator(pos))
            h += interaction(params.predCaptureAndReproduce / 5.0, Predator(left(params.GRIDSIZE)), Predator(pos))
            h += interaction(params.predCaptureAndReproduce / 5.0, Predator(up(params.GRIDSIZE)), Predator(pos))
            h += interaction(params.predCaptureAndReproduce / 5.0, Predator(down(params.GRIDSIZE)), Predator(pos))
            h += interaction(params.predCaptureAndReproduce / 5.0, Predator(pos), Predator(pos))
        }
    }

    override fun toString() = "r($pos)"

    override fun hashCode() = pos*2

    override fun equals(other: Any?) = ((other is Prey) && (pos == other.pos))
}