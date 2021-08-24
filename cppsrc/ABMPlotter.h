//
// Created by daniel on 24/08/2021.
//

#ifndef GLPKTEST_ABMPLOTTER_H
#define GLPKTEST_ABMPLOTTER_H

template<typename AGENT>
class ABMPlotter: Gnuplot {
public:

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

        for (int x = 0; x < PredPreyAgent::GRIDSIZE; ++x) {
            std::vector <HeatRecord> &record = heatData.emplace_back();
            for (int y = 0; y < PredPreyAgent::GRIDSIZE; ++y) {
                double lPrey = analysis[PredPreyAgent(x, y, PredPreyAgent::PREY)];
                double lPred = analysis[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)];
                record.emplace_back(x, y, std::min(lPrey, 1.0) * 200.0, 0.0, std::min(lPred, 1.0) * 200.0);
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
};


#endif //GLPKTEST_ABMPLOTTER_H
