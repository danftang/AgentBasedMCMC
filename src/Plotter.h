//
// Created by daniel on 24/08/2021.
//

#ifndef GLPKTEST_PLOTTER_H
#define GLPKTEST_PLOTTER_H

#include <thread>
#include "gnuplot-iostream/gnuplot-iostream.h"
#include "diagnostics/MultiChainStats.h"

// template specialization to allow valarrays to be recognized


class Plotter: public Gnuplot {
public:

    template<class AGENT, int GRIDSIZE = AGENT::gridsize>
    void plot(const ModelState<AGENT> &realState, const std::valarray<double> &analysis, std::string title="Model state") {
        typedef std::tuple<double, double, double, double, double> HeatRecord;
        std::vector <std::vector<HeatRecord>> heatData;
        std::vector <std::tuple<double, double, double>> pointData;

        // get range of probabilities
        double maxP = -DBL_MAX;
        double minP = DBL_MAX;
//        for(double p : analysis) {
//            if(p > 0.0) {
//                maxP = std::max(maxP, p);
//                minP = std::min(minP, p);
//            }
//        }

        for (int x = 0; x < GRIDSIZE; ++x) {
            for (int y = 0; y < GRIDSIZE; ++y) {
                double gridSquareOccupation = 0.0;
                for(int typeId=0; typeId < AGENT::typeDomainSize; ++typeId) {
                    gridSquareOccupation += analysis[AGENT(x, y, typename AGENT::Type(typeId))];
                }
                gridSquareOccupation /= AGENT::typeDomainSize;
                maxP = std::max(maxP, gridSquareOccupation);
                minP = std::min(minP, gridSquareOccupation);
            }
        }


        std::cout << "Mean ModelState gridsquare occupation range = " << minP << ":" << maxP << std::endl;


        // make realState data
        // map from [0...1] to 0xRRGGBB
        std::function<int(double,double)> colourMap = [](double hue, double saturation) {
            int r = 128.0 * saturation * (1.0 + cos(3.141 * hue));
            int g = 0.0;//255.0 * saturation * sin(3.141 * hue);
            int b = 128.0 * saturation * (1.0 - cos(3.141 * hue));
            return 0x10000*r + 0x100*g + b;
        };
        for (int x = 0; x < GRIDSIZE; ++x) {
            for (int y = 0; y < GRIDSIZE; ++y) {
                double gridsquareMeanOccupation = 0.0;
                int meanType = 0.0;
                for(int typeId=0; typeId < AGENT::typeDomainSize; ++typeId) {
                    int typeOccupation = realState[AGENT(x, y, typename AGENT::Type(typeId))];
                    gridsquareMeanOccupation += typeOccupation;
                    meanType += typeId * typeOccupation;
                }
                if (gridsquareMeanOccupation != 0.0) {
                    double hue = meanType/(AGENT::typeDomainSize*gridsquareMeanOccupation);
                    gridsquareMeanOccupation = std::min(1.0, gridsquareMeanOccupation/AGENT::typeDomainSize);
                    std::cout << gridsquareMeanOccupation << std::endl;
                    pointData.emplace_back(x, y, colourMap(hue,gridsquareMeanOccupation));
                }
            }
        }


        for (int x = 0; x < GRIDSIZE; ++x) {
            std::vector <HeatRecord> &record = heatData.emplace_back();
            for (int y = 0; y < GRIDSIZE; ++y) {
                double gridsquareMeanOccupation = 0.0;
                double meanType = 0.0;
                for(int typeId=0; typeId < AGENT::typeDomainSize; ++typeId) {
                    double typeMeanOccupation = analysis[AGENT(x, y, typename AGENT::Type(typeId))];
                    gridsquareMeanOccupation += typeMeanOccupation;
                    meanType += typeId * typeMeanOccupation;
                }
//                double lPrey = analysis[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREY)];
//                double lPred = analysis[PredPreyAgent<GRIDSIZE>(x, y, PredPreyAgent<GRIDSIZE>::PREDATOR)];
//                lPrey = std::max(lPrey,minP);
//                lPred = std::max(lPred,minP);
//                record.emplace_back(x, y, lPrey * 240.0/maxP, 0.0, lPred * 240.0/maxP);
//                record.emplace_back(x, y, (log(lPrey) - log(minP))* 240.0/(log(maxP)-log(minP)), 0.0, (log(lPred) - log(minP)) * 240.0/(log(maxP)-log(minP)));
                double hue = meanType/(AGENT::typeDomainSize*gridsquareMeanOccupation);
                gridsquareMeanOccupation = std::min(1.0, gridsquareMeanOccupation/(AGENT::typeDomainSize*maxP));
                int rgbColour = colourMap(hue,gridsquareMeanOccupation);
                record.template emplace_back(x,y,rgbColour >> 16, (rgbColour & 0xff00) >> 8, rgbColour & 0xff);
            }
        }

        *this << "set title '" << title << "'\n";
//        *this << "set linetype 1 lc 'red'\n";
//        *this << "set linetype 2 lc 'blue'\n";
//        *this << "set linetype 3 lc 'magenta'\n";

//        *this << "set palette rgbformulae 7,5,15\n";
        *this << "plot [-0.5:" << GRIDSIZE - 0.5 << "][-0.5:" << GRIDSIZE - 0.5 << "] ";
//        *this << "'-' with image notitle, ";
        *this << "'-' with rgbimage notitle";
        if(!pointData.empty()) *this << ", '-' with points pointtype 5 pointsize 0.75 lc rgbcolor variable notitle";
        *this << "\n";
//        *this << "'-' with points pointtype 5 pointsize 0.75 lc palette notitle\n";
        send2d(heatData);
        if(!pointData.empty()) send1d(pointData);
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

    template<int GRIDSIZE, int NTIMESTEPS>
    void animate(Trajectory<PredPreyAgent<GRIDSIZE>,NTIMESTEPS> &predPreyTrajectory, double framesPerSecond) {

        long frameDelay = 1000 / framesPerSecond;
        for(int t=0; t <= NTIMESTEPS; ++t) {
            auto clockTime = std::chrono::steady_clock::now();
            plot(ModelState<PredPreyAgent<GRIDSIZE>>(predPreyTrajectory,t));
            std::this_thread::sleep_until(clockTime + std::chrono::milliseconds(frameDelay));
        }
    }

    template<class AGENT>
    void plot(const MultiChainStats &stats, const ModelState<AGENT> &realEndState = ModelState<AGENT>::zero, bool waitForKeypressToExit = false) {
        // Print MCMC stats
        std::cout << stats;
        // plot autocorrelations
        std::valarray<std::valarray<double>> autocorrelation = stats.autocorrelation();
        int nDimensions = autocorrelation[0].size();
        double xStride = stats.front().varioStride;
        (*this) << "plot ";
        for(int d=1; d<=nDimensions; ++d) (*this) << "'-' using (" << xStride <<  "*$0):" << d << " with lines title 'Statistic " << d << "', ";
        (*this) << "0 with lines notitle\n";
        for(int d=1; d<=nDimensions; ++d) (*this).send1d(autocorrelation);


        // Print scale reduction and effective samples
        std::valarray<double> neff = stats.effectiveSamples();
        std::valarray<double> ineff = (stats.nSamplesPerChain() * 1.0) / neff;
//        auto execTimePerSample = stats.totalSampleTime() * 2.0 / (stats.front().samplerStats.totalProposals() * stats.size());
        auto execTimePerFeasibleSample = stats.totalSampleTime() * 1.0 / (stats.nSamplesPerChain() * stats.nChains());
        std::cout << "Total exec time: " << stats.totalSampleTime() << std::endl;
        std::cout << "Potential scale reduction: " << stats.potentialScaleReduction() << std::endl;
        std::cout << "Actual number of samples per chain: " << stats.nSamplesPerChain() << std::endl;
        std::cout << "Number of chains: " << stats.size() << std::endl;
        std::cout << "Effective number of samples (per chain): " << neff << std::endl;
        std::cout << "Sample inefficiency factor: " << ineff << std::endl << std::endl;
//        std::cout << "Execution time per sample (all threads): " << execTimePerSample << std::endl;
        std::cout << "Execution time per feasible sample (per thread): " << execTimePerFeasibleSample << std::endl;
        std::cout << "Execution time per (worst case) effective sample: " << stats.totalSampleTime()/(neff.min()*stats.nChains()) << std::endl;
        // plot end state

        Plotter endStatePlotter;
        endStatePlotter.plot(realEndState, stats.meanEndState(),"");

        if(waitForKeypressToExit) {
            std::cout << "Press Enter to exit" << std::endl;
            std::cin.get();
        }

    }
};


#endif //GLPKTEST_PLOTTER_H
