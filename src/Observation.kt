import org.apache.commons.math3.fraction.Fraction

// represents a single observation, at a given time
// the constraint is in terms of agent state variables
class Observation(val time: Int, val obs: Constraint<Fraction>) {
}