//
// Created by daniel on 24/08/2021.
//

#ifndef GLPKTEST_PLOTTER_H
#define GLPKTEST_PLOTTER_H

#include "gnuplot-iostream/gnuplot-iostream.h"

// template specialization to allow valarrays to be recognized
template <typename T>
class gnuplotio::ArrayTraits<std::valarray<T>> : public gnuplotio::ArrayTraitsDefaults<T> {
public:
    typedef IteratorRange<const T*, T> range_type;

    static range_type get_range(const std::valarray<T> &arg) {
        return range_type(begin(arg), end(arg));
    }
};


class Plotter: public Gnuplot {
public:

    template<typename AGENT>
    void plot(const ModelState<AGENT> &realState, const std::vector<double> &analysis) {
        typedef std::tuple<double, double, double, double, double> HeatRecord;
        std::vector <std::vector<HeatRecord>> heatData;
        std::vector <std::tuple<double, double, double>> pointData;

//        ModelState <PredPreyAgent> realState = window.realTrajectory.endState();

        for (int x = 0; x < PredPreyAgent::GRIDSIZE; ++x) {
            for (int y = 0; y < PredPreyAgent::GRIDSIZE; ++y) {
                int colour = 2 * (realState[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)] > 0.0)
                             + (realState[PredPreyAgent(x, y, PredPreyAgent::PREY)] > 0.0);
                if (colour != 0)
                    pointData.emplace_back(x, y, colour);
            }
        }

        double maxP = 0.5;
        for (int x = 0; x < PredPreyAgent::GRIDSIZE; ++x) {
            std::vector <HeatRecord> &record = heatData.emplace_back();
            for (int y = 0; y < PredPreyAgent::GRIDSIZE; ++y) {
                double lPrey = analysis[PredPreyAgent(x, y, PredPreyAgent::PREY)];
                double lPred = analysis[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)];
                record.emplace_back(x, y, std::min(lPrey, maxP) * 240.0/maxP, 0.0, std::min(lPred, maxP) * 240.0/maxP);
            }
        }

        *this << "set linetype 1 lc 'red'\n";
        *this << "set linetype 2 lc 'blue'\n";
        *this << "set linetype 3 lc 'magenta'\n";
        *this << "plot [-0.5:" << PredPreyAgent::GRIDSIZE - 0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE - 0.5 << "] ";
        *this << "'-' with rgbimage notitle, ";
        *this << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
        send2d(heatData);
        send1d(pointData);
    }

    template<typename T>
    void heatmap(const T &matrixData) {
//        std::vector<std::vector<double>> stlData(matrixData.size());
//        for(int i=0; i < matrixData.size(); ++i) {
//            stlData[i].reserve(matrixData[i].size());
//            for(int j=0; j < matrixData[i].size(); ++j) {
//                stlData[i].push_back(matrixData[i][j]);
//            }
//        }
//    gp << "plot [-0.5:" << matrixData[0].size() << "][-0.5:" << matrixData.size() << "] ";
        (*this) << "plot '-' matrix with image notitle\n";
        send1d(matrixData);
    }
};


#endif //GLPKTEST_PLOTTER_H
