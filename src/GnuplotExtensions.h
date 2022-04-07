//
// Created by daniel on 21/10/2021.
//

#ifndef GLPKTEST_GNUPLOTEXTENSIONS_H
#define GLPKTEST_GNUPLOTEXTENSIONS_H

//gnuplotio::Gnuplot &operator <<(gnuplotio::Gnuplot &gp, const std::valarray<std::valarray<double>> &twoDData) {
//    std::vector<std::vector<double>> stlData(twoDData.size());
//    for(int i=0; i< twoDData.size(); ++i) {
//        stlData[i].reserve(twoDData[i].size());
//        for(int j=0; j < twoDData[i].size(); ++j) {
//            stlData[i].push_back(twoDData[i][j]);
//        }
//    }
////    gp << "plot [-0.5:" << twoDData[0].size() << "][-0.5:" << twoDData.size() << "] ";
//    gp << "plot '-' matrix with image notitle\n";
//    gp.send1d(stlData);
//    return gp;
//}

#endif //GLPKTEST_GNUPLOTEXTENSIONS_H
