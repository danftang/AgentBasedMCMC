//
// Created by daniel on 24/08/2021.
//

#ifndef GLPKTEST_PLOTTER_H
#define GLPKTEST_PLOTTER_H

#include <thread>
#include "gnuplot-iostream/gnuplot-iostream.h"

// template specialization to allow valarrays to be recognized


class Plotter: public Gnuplot {
public:

    template<int GRIDSIZE>
    void plot(const ModelState<PredPreyAgent<GRIDSIZE>> &realState, const std::valarray<double> &analysis, std::string title="Model state") {
        typedef std::tuple<double, double, double, double, double> HeatRecord;
        std::vector <std::vector<HeatRecord>> heatData;
        std::vector <std::tuple<double, double, double>> pointData;

        for (int x = 0; x < GRIDSIZE; ++x) {
            for (int y = 0; y < GRIDSIZE; ++y) {
                int colour = 2 * (realState[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREDATOR)] > 0.0)
                             + (realState[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREY)] > 0.0);
                if (colour != 0)
                    pointData.emplace_back(x, y, colour);
            }
        }

        double maxP = -DBL_MAX;
        double minP = DBL_MAX;
        for(double p : analysis) {
            if(p > 0.0) {
                maxP = std::max(maxP, p);
                minP = std::min(minP, p);
            }
        }
        std::cout << "P range = " << minP << ":" << maxP << std::endl;

        for (int x = 0; x < GRIDSIZE; ++x) {
            std::vector <HeatRecord> &record = heatData.emplace_back();
            for (int y = 0; y < GRIDSIZE; ++y) {
                double lPrey = analysis[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREY)];
                double lPred = analysis[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREDATOR)];
                lPrey = std::max(lPrey,minP);
                lPred = std::max(lPred,minP);
//                record.emplace_back(x, y, lPrey * 240.0/maxP, 0.0, lPred * 240.0/maxP);
                record.emplace_back(x, y, (log(lPrey) - log(minP))* 240.0/(log(maxP)-log(minP)), 0.0, (log(lPred) - log(minP)) * 240.0/(log(maxP)-log(minP)));
            }
        }

        *this << "set title '" << title << "'\n";
        *this << "set linetype 1 lc 'red'\n";
        *this << "set linetype 2 lc 'blue'\n";
        *this << "set linetype 3 lc 'magenta'\n";
        *this << "plot [-0.5:" << GRIDSIZE - 0.5 << "][-0.5:" << GRIDSIZE - 0.5 << "] ";
        *this << "'-' with rgbimage notitle, ";
        *this << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
        send2d(heatData);
        send1d(pointData);
    }


    template<int GRIDSIZE>
    void plot(const ModelState<PredPreyAgent<GRIDSIZE>> &realState) {
        std::vector <std::tuple<double, double, double>> pointData;

        for (int x = 0; x < GRIDSIZE; ++x) {
            for (int y = 0; y < GRIDSIZE; ++y) {
                int colour = 2 * (realState[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREDATOR)] > 0.0)
                             + (realState[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREY)] > 0.0);
                if (colour != 0)
                    pointData.emplace_back(x, y, colour);
            }
        }

        *this << "set linetype 1 lc 'red'\n";
        *this << "set linetype 2 lc 'blue'\n";
        *this << "set linetype 3 lc 'magenta'\n";
        *this << "plot [-0.5:" << GRIDSIZE - 0.5 << "][-0.5:" << GRIDSIZE - 0.5 << "] ";
        *this << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
        send1d(pointData);
    }


    template<typename T>
    void heatmap(const T &matrixData, double xBegin = 0.0, double xEnd = 0.0, double yBegin = 0.0, double yEnd = 0.0) {
        double xScale = (xEnd == xBegin)?1.0:(xEnd-xBegin)/(matrixData[0].size()-1.0);
        double yScale = (yEnd == yBegin)?1.0:(yEnd-yBegin)/(matrixData.size()-1.0);
        (*this) << "plot '-' using ($1*" <<xScale<< "+" <<xBegin<< "):($2*" <<yScale<< "+"<<yBegin<< "):3 matrix with image notitle\n";
        send1d(matrixData);
    }

    template<int GRIDSIZE>
    void animate(Trajectory<PredPreyAgent<GRIDSIZE>> &predPreyTrajectory, double framesPerSecond) {

        long frameDelay = 1000 / framesPerSecond;
        for(int t=0; t <= predPreyTrajectory.nTimesteps(); ++t) {
            auto clockTime = std::chrono::steady_clock::now();
            plot(ModelState<PredPreyAgent<GRIDSIZE>>(predPreyTrajectory,t));
            std::this_thread::sleep_until(clockTime + std::chrono::milliseconds(frameDelay));
        }
    }
};


#endif //GLPKTEST_PLOTTER_H
