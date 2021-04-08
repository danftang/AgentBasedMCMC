package experiments

import CatAndMouseABM
import GridMapSimplex
import lib.sparseVector.asFractionVector
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import java.util.*

class VertexChooseExpts {

    @Test
    fun chooseVertex() {
        val observations = listOf(CatAndMouseABM.CMObservation(
            CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT),
            1,
            true
        ))
        val constraints = ABMCMC.constraints(CatAndMouseABM, 1, observations)
        println(constraints)
        val simplex = GridMapSimplex(constraints, emptyMap<Int, Fraction>().asFractionVector())
        println(simplex.M)
        println(simplex.X())
        val lexFinal = simplex.getLexFinalBasis()
        println()
        println(simplex.M)
        println(simplex.X())
    }

    fun<T> GridMapSimplex<T>.getLexFinalBasis(): TreeSet<Int> where T : Comparable<T>, T: Number {
        val originalBasis = basicColsByRow.copyOf()
        for(j in nVariables-1 downTo 0) {
            if(!isBasicColumn(j)) {
                val lexEarlierRows = pivotableRows(j, false).filter { basicColsByRow[it] < j }
                lexEarlierRows.firstOrNull()?.also { i ->
                    pivot(i,j)
                }
            }
        }
        initialPivot(originalBasis.asList())
        return TreeSet(basicColsByRow.asList())
    }

    fun<T> GridMapSimplex<T>.gotoLexFirstBasis() where T : Comparable<T>, T: Number {
        for(j in 0 until nVariables) {
            if(!isBasicColumn(j)) {
                val lexEarlierRows = pivotableRows(j, false).filter { basicColsByRow[it] > j }
                lexEarlierRows.firstOrNull()?.also { i ->
                    pivot(i,j)
                }
            }
        }
    }


    // Calculate the probability of choosing the current degeneracy basis, B,
    // given the Lex final basis, F (which is calculated incrementally at each transition).
    //
    // Returns the probability of choosing B
    // if we were to use the following algorithm:
    // Starting with F pivoted-in and B empty, choose a column, v,
    // that
    //  * is higher than the highest variable in B
    //  * is not higher than the lowest variable in F
    //  * has a variable in F as outgoing variable on a pivot.
    // Pivot v in. If there is a choice of outgoing variable in F,
    // choose the lexically earliest. Remove the outgoing
    // variable from F and add v to B. Repeat until B is a complete basis.
    //
    // Not sure about this:
    // Given a tableaux that has B already pivoted-in
    // we can calculate the probability without actually having to do the pivots.
    // Starting with the lowest basic var, v, in B and working up, the options
    // for the next basis are the cols that:
    //  * are higher than the basic var preceeding v in B
    //  * are not greater than the lowest remaining var in F
    //  * have a variable in F that can be a leaving variable on a pivot
    // Since we're dealing with a degenerate basis (i.e. B = 0) if B
    // has a non-zero value in any row that contains a basic variable
    // not before v then, independently of F, it can be pivoted-in without
    // pivoting-out any variables preceeding v in B.
    //
    // To identify the leaving variable in F without pivoting it in, note
    // that the earliest leaving variable in F is the earliest variable in
    // F in the current tableaux that can be pivoted-in with v as leaving variable,
    // which is just the earliest col in F with a non-zero value in v's row.
    // (not true. not sure this is an invariant of the lower pivot state.
    // Instead, it's the earliest member of F that can be pivoted in without displacing
    // a fixed member of B)
    fun<T> GridMapSimplex<T>.calcBasisProbability(lexFinalBasis: TreeSet<Int>): Double where T : Comparable<T>, T: Number {
        TODO()
    }

    // Incremental calculation of final basis and degeneracy probability.
    // Suppose we have already calculated the current final basis and the degeneracy probability
    // What is the effect of adding/removing a row/column?
    // Suppose the rows are ordered by basic variable. Each row is associated with a final column
    // beyond which there is no basis completion, a column in F which is the earliest column
    // that this row's basic var can replace in the hybrid basis and a count of the number
    // of columns that have non-zero values in the submatrix not above this row and between
    // the row's basic var and the final column.
    // The probability is the reciprocal of the product of each row's count.
    // On a pivot, we add and remove rows from the degenerate set, remove a column, add a column
    // and modify some columns.

}