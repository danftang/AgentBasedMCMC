//
// Created by daniel on 06/05/2021.
//

#include <vector>
#include "Experiments.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "SimplexMCMC.h"
#include "agents/PredPreyAgent.h"
#include "PoissonState.h"
#include "DataAssimilation.h"
#include "debug.h"


std::vector<double> Experiments::informationIncrease(int argc, char *argv[]) {
    if(argc != 8) {
        std::cout << "Wrong number of arguments. Should be <GRIDSIZE> <nTimestepsPerWindow> <nWindows> <pPredator> <pPrey> <pMakeObservation> <nSamplesPerWindow> <nBurnInSamples>" << std::endl;
        return std::vector<double>();
    }
    glp_term_out(GLP_OFF); // turn off GLPK terminal output
    int gridsize = atoi(argv[1]);
    int windowSize = atoi(argv[2]);
    int nWindows = atoi(argv[3]);
    double pPredator = atof(argv[4]);         // Poisson prob of predator in each gridsquare at t=0
    double pPrey = atof(argv[5]);             // Poisson prob of prey in each gridsquare at t=0
    double pMakeObservation = atof(argv[6]);  // prob of making an observation of each gridsquare at each timestep
//    double pObserveIfPresent = atof(argv[7]); // prob of detecting an agent given that it is present
    int nSamplesPerWindow = atoi(argv[7]);
    int nBurnInSamples = atoi(argv[8]);
    return informationIncrease(gridsize, windowSize, nWindows, pPredator, pPrey,
                               pMakeObservation, nSamplesPerWindow, nBurnInSamples);
}


std::vector<double> Experiments::informationIncrease(
        int gridsize,
        int windowSize,
        int nWindows,
        double pPredator,         // Poisson prob of predator in each gridsquare at t=0
        double pPrey,             // Poisson prob of prey in each gridsquare at t=0
        double pMakeObservation,  // prob of making an observation of each gridsquare at each timestep
//        double pObserveIfPresent, // prob of detecting an agent given that it is present
        int nSamplesPerWindow,
        int nBurnInSamples
        ) {
    glp_term_out(GLP_OFF); // turn off GLPK terminal output
    PredPreyAgent::GRIDSIZE = gridsize;

    PoissonState<PredPreyAgent> startState([&](const PredPreyAgent &agent) {
        return agent.type()==PredPreyAgent::PREDATOR?pPredator:pPrey;
//        return agent.type() == PredPreyAgent::PREDATOR ?
//               (pPredator *(agent.xPosition() < PredPreyAgent::GRIDSIZE / 2)) :
//               (pPrey * (agent.xPosition() >= PredPreyAgent::GRIDSIZE / 2));
    });

    auto observationOperator = [=](const Trajectory<PredPreyAgent> &trajectory) {
        return Observation<PredPreyAgent>::generateObservations(trajectory, pMakeObservation);
    };

    DataAssimilation<PredPreyAgent> assimilation(nWindows, windowSize, startState, observationOperator,
                                                 nSamplesPerWindow, nBurnInSamples);
    return assimilation.calculateInformationGain();
}


void Experiments::PredPreyAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 8;
    constexpr int windowSize = 2;
    constexpr int nWindows = 1;
    constexpr double pPredator = 0.16;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 0.32;//2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.05;//0.04;    // prob of making an observation of each gridsquare at each timestep
//    constexpr double pObserveIfPresent = 0.999; // 0.9;
    constexpr int nSamplesPerWindow = 25000; //250000;
    constexpr int nBurninSamples = 1000;
//    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////
    PoissonState<PredPreyAgent> prior([](const PredPreyAgent &agent) {
        return agent.type()==PredPreyAgent::PREDATOR?pPredator:pPrey;
//        return agent.type()==PredPreyAgent::PREDATOR?(pPredator*(agent.xPosition() < PredPreyAgent::GRIDSIZE/2)):(pPrey*(agent.xPosition() >= PredPreyAgent::GRIDSIZE/2));
    });

    auto observationOperator = [=](const Trajectory<PredPreyAgent> &trajectory) {
        return Observation<PredPreyAgent>::generateObservations(trajectory, pMakeObservation);
    };

    DataAssimilation<PredPreyAgent> assimilation(prior, observationOperator);


    for(int w=0; w<nWindows; ++w) {
        assimilation.addWindow(windowSize, nSamplesPerWindow, nBurninSamples);
        Gnuplot gp;
        gp << assimilation.windows[w];
    }

//    std::vector<PoissonState<PredPreyAgent>> priors = assimilation.priorWindows(1000);
//    for(int w=0; w<nWindows; ++w) {
//        Gnuplot gp;
//        gp << priors[w];
//    }
}




void Experiments::PredPreyExpt() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 8;
    constexpr int nTimesteps = 4;
    constexpr double pPredator = 0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.02;    // prob of making an observation of each gridsquare at each timestep
//    constexpr double pObserveIfPresent = 0.95; // 0.9;
    constexpr int nSamples = 10000; //250000;
    constexpr int nBurnInSamples = 1000;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////
    PoissonState<PredPreyAgent> prior([](const PredPreyAgent &agent) {
        return agent.type()==PredPreyAgent::PREDATOR?pPredator:pPrey;
    });

    ModelState<PredPreyAgent> startState = prior.nextSample();

    std::cout << "Start state: " << startState << std::endl;
    Trajectory<PredPreyAgent> realTrajectory(nTimesteps, startState);
    std::vector<Observation<PredPreyAgent>> observations =
            Observation<PredPreyAgent>::generateObservations(realTrajectory, pMakeObservation);
    std::cout << "Real trajectory: " << glp::SparseVec(realTrajectory) << std::endl;
    std::cout << "Observations: " << observations << std::endl;

    AssimilationWindow<PredPreyAgent> window(realTrajectory, prior, observations, nSamples, nBurnInSamples);


    ////////////////////////////////////////// SHOW RESULTS ///////////////////////////////////////////
    std::cout << std::endl;
    Gnuplot gp;
    plotHeatMap(gp, window.analysis, window.realTrajectory.endState());
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Experiments::CatMouseExpt() {
    CatMouseAgent leftCat(CatMouseAgent::Type::CAT, CatMouseAgent::Position::LEFT);
    auto observations = std::vector({ Observation(1,leftCat, 1) });

    ABMProblem<CatMouseAgent> abm(2, observations, [](const Trajectory<CatMouseAgent> &trajectory) {
        return 0.0;
    });

    std::cout << abm << std::endl;

    // calculate initial trajectory
    abm.cpxBasis();
    abm.simplex();
    std::cout << "LP relaxation initial trajectory: " << abm.primalSolution() << std::endl;
    abm.intOpt();
    std::cout << "MIP initial trajectory: " << abm.mipSolution() << std::endl;
    abm.warmUp();

//    Trajectory<CatMouseAgent> initialTrajectory;
//    initialTrajectory.add(Event(0,leftCat, CatMouseAgent::STAYPUT),1.0);
//    std::cout << "Initial trajectory is" << std::endl;
//    std::cout << initialTrajectory << std::endl;
//    abm.stdBasis();
//    abm.warmUp();
    abm.advBasis();
    abm.warmUp();
    SimplexMCMC mcmc(abm, abm.logProbFunc());
    std::cout << mcmc.X() << std::endl;
    for(int n=0; n<100; ++n) {
        mcmc.nextSample();
//        mcmc.randomWalk();
        std::cout << n << "  " << mcmc.X() << std::endl;
        std::cout << "Valid = " << abm.isValidSolution(mcmc.X()) << std::endl;
        assert(abm.isValidSolution(mcmc.X()));
    }

}


void Experiments::RandomWalk() {
    using glp::X;
    glp::Problem myProb;

    myProb.addConstraint(1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 100.0);
    myProb.addConstraint(10.0*X(1) + 4.0*X(2) + 5.0*X(3) <= 600.0);
    myProb.addConstraint(2.0*X(1) + 2.0*X(2) + 6.0*X(3) <= 300.0);
    myProb.addConstraint(0.0 <= 1.0*X(1));
    myProb.addConstraint(0.0 <= 1.0*X(2));
    myProb.addConstraint(0.0 <= 1.0*X(3));
    myProb.setObjective(10.0*X(1) + 6.0*X(2) + 4.0*X(3));
    myProb.stdBasis();
    myProb.warmUp();
    std::cout << myProb;

    SimplexMCMC myMCMC(myProb, nullPMF);

    std::cout << myMCMC << std::endl;

    for(int sample=0; sample < 5; ++sample) {
        myMCMC.randomWalk();
        std::cout << myMCMC << std::endl;
        std::cout << "Sample is: " << myMCMC.X() << std::endl;
    }

//    mySimplex.pivot(3,3);
//    std::cout << mySimplex << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//void Experiments::GnuplotTest() {
//    Gnuplot gp; //(stdout);
//    auto state = ModelState<PredPreyAgent>::randomPoissonState([](const PredPreyAgent &agent) {
//        if(agent.type() == PredPreyAgent::PREDATOR) return 0.08;
//        return 0.16;
//    });
//
//    plotHeatMap(gp, state, state);
//}

Gnuplot &Experiments::plotHeatMap(Gnuplot &gp, const PoissonState<PredPreyAgent> &aggregateState,
                                  const ModelState<PredPreyAgent> &realState) {
    typedef std::tuple<double,double,double,double,double> HeatRecord;
    std::vector<std::vector<HeatRecord>> heatData;
    std::vector<std::tuple<double,double,double>> pointData;

//    for(auto [agent, occupancy] : realState) {
//        pointData.emplace_back(agent.xPosition(), agent.yPosition(), agent.type()==PredPreyAgent::PREY?1:2);
//    }

    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
            int colour = 2*(realState[PredPreyAgent(x,y,PredPreyAgent::PREDATOR)]>0.0)
                    + (realState[PredPreyAgent(x,y,PredPreyAgent::PREY)]>0.0);
            if(colour != 0)
                pointData.emplace_back(x, y, colour);
        }
    }


//    double predMaxOccupancy = 0.0;
//    double preyMaxOccupancy = 0.0;
//    for(auto [agent, occupancy] : aggregateState) {
//        if(agent.type() == PredPreyAgent::PREDATOR) {
//            if(occupancy > predMaxOccupancy) predMaxOccupancy = occupancy;
//        } else {
//            if(occupancy > preyMaxOccupancy) preyMaxOccupancy = occupancy;
//        }
//    }
    double maxOccupancy = 0.0;
    for(double occupancy : aggregateState.stateCounts) {
        if(occupancy > maxOccupancy) maxOccupancy = occupancy;
    }
//    std::cout << "Max occupancy = " << maxOccupancy << " at " << maxState << std::endl;

//    double predScale = 200.0/log(predMaxOccupancy + 1.0);
//    double preyScale = 200.0/log(preyMaxOccupancy + 1.0);
    double predScale = 200.0/maxOccupancy; // predMaxOccupancy;
    double preyScale = 200.0/maxOccupancy; // preyMaxOccupancy;
    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
        std::vector<HeatRecord> &record = heatData.emplace_back();
        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
            double nPrey = aggregateState.stateCounts[PredPreyAgent(x, y, PredPreyAgent::PREY)];
            double nPred = aggregateState.stateCounts[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)];
//            record.emplace_back(x, y, log(nPrey + 1.0) * preyScale, 0.0, log(nPred + 1.0) * predScale);
            record.emplace_back(x,y,nPrey*preyScale,0.0,nPred*predScale);
        }
    }

    gp << "set linetype 1 lc 'red'\n";
    gp << "set linetype 2 lc 'blue'\n";
    gp << "set linetype 3 lc 'magenta'\n";
    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
    gp << "'-' with rgbimage notitle, ";
    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
    gp.send2d(heatData);
    gp.send1d(pointData);
    return gp;
}


Gnuplot &Experiments::plotAgents(Gnuplot &gp, const ModelState<PredPreyAgent> &state) {
    std::vector<std::tuple<double,double,double>> pointData;
    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
            int colour = 2*(state[PredPreyAgent(x,y,PredPreyAgent::PREDATOR)]>0.0)
                         + (state[PredPreyAgent(x,y,PredPreyAgent::PREY)]>0.0);
            if(colour != 0)
                pointData.emplace_back(x, y, colour);
        }
    }

    gp << "set linetype 1 lc 'red'\n";
    gp << "set linetype 2 lc 'blue'\n";
    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable\n";
    gp.send1d(pointData);
    return gp;
}

