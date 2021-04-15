//#include <iostream>
#include <iostream>
#include "glpk.h"
#include "spxprob.h"
#include "GlpProblem.h"
#include "GlpTableau.h"
#include "GlpSimplex.h"


int main() {
    glp_prob *lp;
    int ia[1+1000], ja[1+1000];
    double ar[1+1000], z, x1, x2, x3;
    s1:   lp = glp_create_prob();
    s2:   glp_set_prob_name(lp, "sample");
    s3:   glp_set_obj_dir(lp, GLP_MIN);
    s4:   glp_add_rows(lp, 3);
    s5:   glp_set_row_name(lp, 1, "p");
    s6:   glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 100.0);
    s7:   glp_set_row_name(lp, 2, "q");
    s8:   glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 600.0);
    s9:   glp_set_row_name(lp, 3, "r");
    s10:  glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 300.0);
    s11:  glp_add_cols(lp, 3);
    s12:  glp_set_col_name(lp, 1, "x1");
    s13:  glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);
    s14:  glp_set_obj_coef(lp, 1, 10.0);
    s15:  glp_set_col_name(lp, 2, "x2");
    s16:  glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);
    s17:  glp_set_obj_coef(lp, 2, 6.0);
    s18:  glp_set_col_name(lp, 3, "x3");
    s19:  glp_set_col_bnds(lp, 3, GLP_LO, 0.0, 0.0);
    s20:  glp_set_obj_coef(lp, 3, 4.0);
    s21:  ia[1] = 1, ja[1] = 1, ar[1] =  1.0; /* a[1,1] =  1 */
    s22:  ia[2] = 1, ja[2] = 2, ar[2] =  1.0; /* a[1,2] =  1 */
    s23:  ia[3] = 1, ja[3] = 3, ar[3] =  1.0; /* a[1,3] =  1 */
    s24:  ia[4] = 2, ja[4] = 1, ar[4] = 10.0; /* a[2,1] = 10 */
    s26:  ia[6] = 2, ja[6] = 2, ar[6] =  4.0; /* a[2,2] =  4 */
    s28:  ia[8] = 2, ja[8] = 3, ar[8] =  5.0; /* a[2,3] =  5 */
    s25:  ia[5] = 3, ja[5] = 1, ar[5] =  2.0; /* a[3,1] =  2 */
    s27:  ia[7] = 3, ja[7] = 2, ar[7] =  2.0; /* a[3,2] =  2 */
    s29:  ia[9] = 3, ja[9] = 3, ar[9] =  1.0; /* a[3,3] =  6 */
    s30:  glp_load_matrix(lp, 9, ia, ja, ar);
    glp_std_basis(lp);
    glp_warm_up(lp);
//    glp_factorize(lp);
//    s31:  glp_simplex(lp, NULL);
//    s32:  z = glp_get_obj_val(lp);
//    s33:  x1 = glp_get_col_prim(lp, 1);
//    s34:  x2 = glp_get_col_prim(lp, 2);
//    s35:  x3 = glp_get_col_prim(lp, 3);
//    s36:  printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n",
//                 z, x1, x2, x3);


    GlpProblem myProb(lp);
//    GlpTableau myTableau(lp);
    GlpSimplex mySimplex(lp);
//    SparseVec c(myTableau.nRows());

    std::cout << myProb;
//    std::cout << myTableau;
    std::cout << mySimplex << std::endl;

    mySimplex.pivot(3,3,true);

    std::cout << mySimplex << std::endl;

//    std::cout << "Columns" << std::endl;
//
//    for(int k=0; k<6; ++k) {
//        myTableau.col(k,c);
//        std::cout << c << std::endl;
//    }
//
////    for(int i=0; i<3; ++i) {
////        std::cout << glp_eval_tab_col(lp,k,) << std::endl;
////    }
//
//
////    glp_set_col_stat(lp, 3, GLP_BS);
////    glp_set_row_stat(lp, 3, GLP_NU);
//    myTableau.pivot(2,5);
////    glp_unscale_prob(lp);
////
////    std::cout << std::endl;
////    std::cout << myTableau;
//
////    lp->valid = 0;
////    glp_factorize(lp);
////    std::cout << std::endl;
//////    glp_unscale_prob(lp);
//    std::cout << std::endl;
//    std::cout << myTableau;
////
//    for(int k=0; k<6; ++k) {
//        myTableau.col(k,c);
//        std::cout << c << std::endl;
//    }
////
////


s37:  glp_delete_prob(lp);

//    glp_print_tableau(lp);
//    bfd_update(lp->bfd, 1, )


    return 0;
}


//void buildSimplex(glp_prob *P, const glp_smcp *parm) {
//    SPXLP lp;
//    int *map;
//    spx_init_lp(&lp, P, parm->excl);
//    spx_alloc_lp(&lp);
//    map = talloc(1+P->m+P->n, int);
//    spx_build_lp(&lp, P, parm->excl, parm->shift, map);
//    spx_build_basis(&lp, P, map);
//
//}
//
//void getBf(glp_prob *P, const glp_smcp *parm) {
//    BFD *bfd = P->bfd;
//    int idx[2];
//    double vals[2];
//    bfd_update(bfd, 1, 2, idx, vals);
//    glp_factorize(P);
//
//    FVS sparseVec;
//    sparseVec.n = 1000;
//    sparseVec.nnz = 1;
//    sparseVec.ind = idx;
//    sparseVec.vec = vals;
//    bfd_ftran_s(bfd, &sparseVec);
//}
