import com.google.ortools.linearsolver.MPConstraint
import com.google.ortools.linearsolver.MPSolver
import com.google.ortools.linearsolver.MPVariable
import com.google.ortools.sat.CpModel
import com.google.ortools.sat.IntVar
import lib.Multiset
import lib.MutableMultiset

interface UnknownModelState<AGENT> {
    val potentialConsequencesFootprint: Set<AGENT>
    val committedConsequences: Multiset<AGENT>

    fun getEndStateVriable(model: CpModel, agent: AGENT): IntVar
    fun getEndStateBinaryVriable(model: CpModel, agent: AGENT): IntVar
}
