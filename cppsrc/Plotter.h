//
// Created by daniel on 24/08/2021.
//

#ifndef GLPKTEST_PLOTTER_H
#define GLPKTEST_PLOTTER_H

#include "gnuplot-iostream/gnuplot-iostream.h"

// template specialization to allow valarrays to be recognized


class Plotter: public Gnuplot {
public:

    template<int GRIDSIZE>
    void plot(const ModelState<PredPreyAgent<GRIDSIZE>> &realState, const std::vector<double> &analysis) {
        typedef std::tuple<double, double, double, double, double> HeatRecord;
        std::vector <std::vector<HeatRecord>> heatData;
        std::vector <std::tuple<double, double, double>> pointData;

//        ModelState <PredPreyAgent> realState = window.realTrajectory.endState();

        for (int x = 0; x < GRIDSIZE; ++x) {
            for (int y = 0; y < GRIDSIZE; ++y) {
                int colour = 2 * (realState[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREDATOR)] > 0.0)
                             + (realState[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREY)] > 0.0);
                if (colour != 0)
                    pointData.emplace_back(x, y, colour);
            }
        }

        double maxP = 0.5;
        for (int x = 0; x < GRIDSIZE; ++x) {
            std::vector <HeatRecord> &record = heatData.emplace_back();
            for (int y = 0; y < GRIDSIZE; ++y) {
                double lPrey = analysis[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREY)];
                double lPred = analysis[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREDATOR)];
                record.emplace_back(x, y, std::min(lPrey, maxP) * 240.0/maxP, 0.0, std::min(lPred, maxP) * 240.0/maxP);
            }
        }

        *this << "set linetype 1 lc 'red'\n";
        *this << "set linetype 2 lc 'blue'\n";
        *this << "set linetype 3 lc 'magenta'\n";
        *this << "plot [-0.5:" << GRIDSIZE - 0.5 << "][-0.5:" << GRIDSIZE - 0.5 << "] ";
        *this << "'-' with rgbimage notitle, ";
        *this << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
        send2d(heatData);
        send1d(pointData);
    }

    template<typename T>
    void heatmap(const T &matrixData, double xBegin = 0.0, double xEnd = 0.0, double yBegin = 0.0, double yEnd = 0.0) {
        double xScale = (xEnd == xBegin)?1.0:(xEnd-xBegin)/(matrixData[0].size()-1.0);
        double yScale = (yEnd == yBegin)?1.0:(yEnd-yBegin)/(matrixData.size()-1.0);
        (*this) << "plot '-' using ($1*" <<xScale<< "+" <<xBegin<< "):($2*" <<yScale<< "+"<<yBegin<< "):3 matrix with image notitle\n";
        send1d(matrixData);
    }
};


#endif //GLPKTEST_PLOTTER_H
