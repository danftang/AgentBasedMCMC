//
// Created by daniel on 06/05/2021.
//

#include <vector>
#include "Experiments.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "SimplexMCMC.h"
#include "agents/PredPreyAgent.h"

void Experiments::PredPreyExpt() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 10;
    constexpr int nTimesteps = 32;
    constexpr double pPredator = 0.04;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 4.0*pPredator;//2.0*pPredator;     // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.02;    // prob of making an observation of each gridsquare at each timestep
    constexpr int nSamples = 1; //250000;
    constexpr int plotTimestep = 0; //nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////
    ModelState<PredPreyAgent> startState = ModelState<PredPreyAgent>::randomPoissonState([](const PredPreyAgent &agent) {
        if(agent.type() == PredPreyAgent::PREDATOR) return pPredator;
        return pPrey;
    });
    std::cout << "Start state: " << startState << std::endl;
    auto [observations, realTrajectory] =
            Observation<PredPreyAgent>::generateObservations(startState, nTimesteps, pMakeObservation);
    ABMProblem<PredPreyAgent> abm(nTimesteps, observations, [](const Trajectory<PredPreyAgent> &trajectory) {
        ModelState<PredPreyAgent> startState = trajectory(0);
        double logP = 0.0;
        for(auto [agent, occupation]: startState) {
            double k = fabs(occupation);
            double lambda = agent.type()==PredPreyAgent::PREDATOR?pPredator:pPrey;
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
//    mcmc.setLPState(realTrajectory);

    // solve by phase1
    std::cout << "Starting phase 1 in state: " << glp::SparseVec(mcmc.X()) << std::endl;
    mcmc.findFeasibleStartPoint();

    ////////////////////////////////////////// DO SANITY CHECKS ////////////////////////////////////////
    // Check initial basis contains no fixed vars and all auxiliaries are in the basis
    for(int k=1; k<=mcmc.nVars(); ++k) {
        if(mcmc.l[k] == mcmc.u[k]) {
            std::cout << "WARNING: Fixed variable in the Simplex" << k << std::endl;
            return;
        }
        if(mcmc.kSimTokProb[k] > abm.nConstraints() && (mcmc.l[k] != 0.0 || mcmc.u[k] != 1.0 )) {
            std::cout << "WARNING: non-binary structural var " << k << std::endl;
            return;
        }
    }
    for(int j=1; j<=mcmc.nNonBasic(); ++j) {
        int k = mcmc.head[mcmc.nBasic()+j];
        if(mcmc.kSimTokProb[k] <= abm.nConstraints()) {
            std::cout << "WARNING: Non-basic auxiliary variable " << j << std::endl;
            return;
        }
        if(mcmc.l[k] == -DBL_MAX || mcmc.u[k] == DBL_MAX) {
            std::cout << "WARNING: Non-basic variable with infinite bound " << j << std::endl;
            return;
        }
    }
    std::cout << "Starting with initial sample:" << std::endl;
    std::cout << glp::SparseVec(mcmc.X()) << std::endl;
    assert(abm.isValidSolution(mcmc.X()));

    ////////////////////////////////////////// DO SAMPLING ///////////////////////////////////////////
    ModelState<PredPreyAgent> meanState;
    for(int n=0; n<nSamples; ++n) {
        mcmc.nextSample();
        if(nSamples<1000 || n%100 == 1) {
            std::cout << "Sample " << n << std::endl;
            assert(abm.isValidSolution(mcmc.X()));
//            std::cout << "Sample " << n << " : " << glp::SparseVec(mcmc.X()) << std::endl;
        }
        meanState += (reinterpret_cast<const Trajectory<PredPreyAgent> *>(&mcmc.X()))->operator()(plotTimestep); // TODO: this is ugly
    }


    ////////////////////////////////////////// SHOW RESULTS ///////////////////////////////////////////
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
    auto observations = std::vector(
            {Observation(
                    1,
                    leftCat,
                    1,
                    0.95)
            });

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
void Experiments::GnuplotTest() {
    Gnuplot gp; //(stdout);
    auto state = ModelState<PredPreyAgent>::randomPoissonState([](const PredPreyAgent &agent) {
        if(agent.type() == PredPreyAgent::PREDATOR) return 0.08;
        return 0.16;
    });
    plotHeatMap(gp, state, state);
}

Gnuplot &Experiments::plotHeatMap(Gnuplot &gp, const ModelState<PredPreyAgent> &aggregateState,
                                  const ModelState<PredPreyAgent> &realState) {
    typedef std::tuple<double,double,double,double,double> HeatRecord;
    std::vector<std::vector<HeatRecord>> heatData;
    std::vector<std::tuple<double,double,double>> pointData;

    for(auto [agent, occupancy] : realState) {
        pointData.emplace_back(agent.xPosition(), agent.yPosition(), agent.type()==PredPreyAgent::PREY?1:2);
    }

    double predMaxOccupancy = 0.0;
    double preyMaxOccupancy = 0.0;
    for(auto [agent, occupancy] : aggregateState) {
        if(agent.type() == PredPreyAgent::PREDATOR) {
            if(occupancy > predMaxOccupancy) predMaxOccupancy = occupancy;
        } else {
            if(occupancy > preyMaxOccupancy) preyMaxOccupancy = occupancy;
        }
    }

//    double predScale = 200.0/log(predMaxOccupancy + 1.0);
//    double preyScale = 200.0/log(preyMaxOccupancy + 1.0);
    double predScale = 224.0/predMaxOccupancy;
    double preyScale = 224.0/preyMaxOccupancy;
    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
        std::vector<HeatRecord> &record = heatData.emplace_back();
        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
            double nPrey = aggregateState[PredPreyAgent(x, y, PredPreyAgent::PREY)];
            double nPred = aggregateState[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)];
//            record.emplace_back(x, y, log(nPrey + 1.0) * preyScale, 0.0, log(nPred + 1.0) * predScale);
            record.emplace_back(x,y,nPrey*preyScale,0.0,nPred*predScale);
        }
    }

    gp << "set linetype 1 lc 'red'\n";
    gp << "set linetype 2 lc 'blue'\n";
    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
    gp << "'-' with rgbimage notitle, ";
    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
    gp.send2d(heatData);
    gp.send1d(pointData);
    return gp;
}


Gnuplot &Experiments::plotAgents(Gnuplot &gp, const ModelState<PredPreyAgent> &state) {
    std::vector<std::tuple<double,double,double>> pointData;
    for(auto [agent, occupation]: state) {
        pointData.emplace_back(agent.xPosition(), agent.yPosition(), agent.type()==PredPreyAgent::PREY?1:2);
    }
    gp << "set linetype 1 lc 'red'\n";
    gp << "set linetype 2 lc 'blue'\n";
    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable\n";
    gp.send1d(pointData);
    return gp;
}
