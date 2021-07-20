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
    std::vector<double> informationGain;
    if(argc != 9) {
        std::cout << "Wrong number of arguments. Should be <GRIDSIZE> <nTimestepsPerWindow> <nWindows> <pPredator> <pPrey> <pMakeObservation> <pObserveIfPresent> <nSamplesPerWindow>" << std::endl;
    } else {
        glp_term_out(GLP_OFF); // turn off GLPK terminal output
        PredPreyAgent::GRIDSIZE = atoi(argv[1]);
        int windowSize = atoi(argv[2]);
        int nWindows = atoi(argv[3]);
        double pPredator = atof(argv[4]);//0.08;          // Poisson prob of predator in each gridsquare at t=0
        double pPrey = atof(argv[5]);//2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
        double pMakeObservation = atof(
                argv[6]);//0.04;    // prob of making an observation of each gridsquare at each timestep
        double pObserveIfPresent = atof(argv[7]); // 0.9;
        int nSamplesPerWindow = atoi(argv[8]); //250000;


        PoissonState<PredPreyAgent> startState([&](const PredPreyAgent &agent) {
            return agent.type() == PredPreyAgent::PREDATOR ? (pPredator *
                                                              (agent.xPosition() < PredPreyAgent::GRIDSIZE / 2)) : (
                           pPrey * (agent.xPosition() >= PredPreyAgent::GRIDSIZE / 2));
        });

        DataAssimilation<PredPreyAgent> assimilation(nWindows, windowSize, startState, pMakeObservation, pObserveIfPresent, nSamplesPerWindow);


//        std::vector<PoissonState<PredPreyAgent>> priors = DataAssimilation<PredPreyAgent>::prior(nWindows, windowSize,
//                                                                                                nSamplesPerWindow,
//                                                                                                startState);
//
//        DataAssimilation<PredPreyAgent>::runAndObserve(nWindows, windowSize, pMakeObservation, pObserveIfPresent, startState.sample());
//
//        std::vector<PoissonState<PredPreyAgent>> posterior = DataAssimilation<PredPreyAgent>::posterior(nWindows,
//                                                                                                        windowSize,
//                                                                                                        nSamplesPerWindow,
//                                                                                                        startState);
//
//        for (int t = 0; t < nWindows; ++t) {
//            informationGain.push_back(DataAssimilation<PredPreyAgent>::informationGain());
//        }
    }
    return informationGain;
}


void Experiments::PredPreyAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 12;
    constexpr int windowSize = 6;
    constexpr int nWindows = 2;
    constexpr double pPredator = 0.16;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 0.32;//2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.25;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 0.95; // 0.9;
    constexpr int nSamplesPerWindow = 50000; //250000;
//    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////
    PoissonState<PredPreyAgent> poissonModelState([](const PredPreyAgent &agent) {
//        return agent.type()==PredPreyAgent::PREDATOR?pPredator:pPrey;
        return agent.type()==PredPreyAgent::PREDATOR?(pPredator*(agent.xPosition() < PredPreyAgent::GRIDSIZE/2)):(pPrey*(agent.xPosition() >= PredPreyAgent::GRIDSIZE/2));
    });

    DataAssimilation<PredPreyAgent> assimilation(nWindows, windowSize, poissonModelState, pMakeObservation, pObserveIfPresent, nSamplesPerWindow);


//    ModelState<PredPreyAgent> realState = poissonModelState.sample();
//
//    Gnuplot gp0;
//    plotHeatMap(gp0, poissonModelState, realState);
//
//    for(int window=0; window<nWindows; ++window) {
//        debug(std::cout << "Real state: " << realState << std::endl);
//        debug(std::cout << "Poission state: " << poissonModelState << std::endl);
//        auto [observations, realTrajectory] =
//        Observation<PredPreyAgent>::generateObservations(realState, windowSize, pMakeObservation, pObserveIfPresent);
////        std::cout << "Real trajectory: " << glp::SparseVec(realTrajectory) << std::endl;
//        debug(std::cout << "Observations: " << observations << std::endl);
//
//        poissonModelState = AssimilationWindow<PredPreyAgent>(
//                windowSize,
//                observations,
//                poissonModelState,
//                nSamplesPerWindow
//        );
//        realState = realTrajectory(windowSize);
//        Gnuplot gp;
//        plotHeatMap(gp, poissonModelState, realState);
//    }

}


void Experiments::PredPreyExpt() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 8;
    constexpr int nTimesteps = 4;
    constexpr double pPredator = 0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.02;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 0.95; // 0.9;
    constexpr int nSamples = 10000; //250000;
    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////
    ModelState<PredPreyAgent> startState = ModelState<PredPreyAgent>::randomPoissonState([](const PredPreyAgent &agent) {
        if(agent.type() == PredPreyAgent::PREDATOR) return pPredator;
        return pPrey;
    });
    std::cout << "Start state: " << startState << std::endl;
    auto [observations, realTrajectory] =
            Observation<PredPreyAgent>::generateObservations(startState, nTimesteps, pMakeObservation, pObserveIfPresent);
    ABMProblem<PredPreyAgent> abm(nTimesteps, observations, [](const Trajectory<PredPreyAgent> &trajectory) {
        ModelState<PredPreyAgent> startState = trajectory(0);
        double logP = 0.0;
        for(int agentId=0; agentId < PredPreyAgent::domainSize(); ++agentId) {
            double k = fabs(startState[agentId]);
            double lambda = PredPreyAgent(agentId).type()==PredPreyAgent::PREDATOR?pPredator:pPrey;
            logP += k*log(lambda) - lambda - lgamma(k+1); // log of Poisson
        }
        return logP;
    });
    std::cout << "Real trajectory: " << glp::SparseVec(realTrajectory) << std::endl;
    std::cout << "Observations: " << observations << std::endl;

    SimplexMCMC mcmc(abm, abm.logProbFunc());

    //////////////////////////////////// FIND INITIAL SOLUTION ////////////////////////////////////////
//    // solve by simplex
//    abm.simplex();
//    std::vector<double> initSol = abm.primalSolution();
//    assert(abm.isValidSolution(initSol));
//    std::cout << "Found initial solution: " << glp::SparseVec(initSol) << std::endl;
//    mcmc.setLPState(initSol);

    // use real trajectory as initial state
    mcmc.setLPState(realTrajectory);

    // solve by phase1
//    std::cout << "Starting phase 1 in state: " << glp::SparseVec(mcmc.X()) << std::endl;
//    mcmc.findFeasibleStartPoint();

    std::cout << "Starting with initial sample:" << std::endl;
    std::cout << glp::SparseVec(mcmc.X()) << std::endl;

    ////////////////////////////////////////// DO SANITY CHECKS ////////////////////////////////////////
    assert(mcmc.abmSanityChecks());


    ////////////////////////////////////////// DO SAMPLING ///////////////////////////////////////////
    PoissonState<PredPreyAgent> meanState;
    for(int n=0; n<nSamples; ++n) {
        mcmc.nextSample();
        if(nSamples<1000 || n%100 == 1) {
            std::cout << "Sample " << n << std::endl;
            assert(abm.isValidSolution(mcmc.X()));
//            std::cout << "Sample " << n << " : " << glp::SparseVec(mcmc.X()) << std::endl;
        }
        const Trajectory<PredPreyAgent> &trajectory = reinterpret_cast<const Trajectory<PredPreyAgent> &>(mcmc.X());
        meanState += trajectory(plotTimestep);
    }


    ////////////////////////////////////////// SHOW RESULTS ///////////////////////////////////////////
    std::cout << std::endl;
    std::cout << "Feasible sample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl;
    std::cout << "Infeasible sample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
    std::cout << "Mean state:\n" << meanState << std::endl;
    std::cout << "Real state:\n" << realTrajectory(plotTimestep) << std::endl;
    std::cout << "Observations: " << observations << std::endl;
    Gnuplot gp;
//    StateTrajectory<PredPreyAgent> realStateTrajectory(realTrajectory);
    plotHeatMap(gp, meanState, realTrajectory(plotTimestep));
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Experiments::CatMouseExpt() {
    CatMouseAgent leftCat(CatMouseAgent::Type::CAT, CatMouseAgent::Position::LEFT);
    auto observations = std::vector({ Observation(State<CatMouseAgent>(1,leftCat), 1, 0.95) });

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

