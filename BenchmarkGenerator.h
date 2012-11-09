// $Id: BenchmarkGenerator.h 19 2011-06-17 08:20:59Z jbao $

#ifndef BENCHMARKGENERATOR_H
#define BENCHMARKGENERATOR_H

#include "nr/nr.h"
#include "GeneNetwork.h"
#include "SteadyStateExperiment.h"
#include "TimeSeriesExperiment.h"

class BenchmarkGenerator {

public:
	BenchmarkGenerator(const GeneNetwork& grn);
	~BenchmarkGenerator();
	void generateGoldStandard();

private:
	/** The gene network */
	GeneNetwork grn_;
	/** The steady-state experiments */
	std::vector<SteadyStateExperiment*> steadyStateExperiments_;
	/** If both GnwSettings.simulateODE_ and simulateSDE_ are set, these are the simulations using ODEs */
	//std::vector<SteadyStateExperiment> steadyStateExperimentsODE_;
	/** The time-series experiments */
	std::vector<TimeSeriesExperiment> timeSeriesExperiments_;
	/** If both GnwSettings.simulateODE_ and simulateSDE_ are set, these are the simulations using ODEs */
	//ArrayList<TimeSeriesExperiment> timeSeriesExperimentsODE_;
	/** The wild-type obtained using ODEs is often used as initial condition */
	Vec_DP wildTypeODE_;
	/** The perturbations applied to the time-series */
	//Perturbation timeSeriesPerturbations_;
	
	void runAll();
	void runSteadyStateExperiments();
	void runTimeSeriesExperiments() throw(std::exception);
	void constructInitialConditionFromWildType(Mat_I_DP& perturb, Vec_O_DP& xy0);
	void createAndRunTsExperiment(Perturbation *perturbation, bool restoreWildTypeAtHalftime, std::string label, Vec_DP xy0);
	void addExperimentalNoise();
	void normalize();
	
};

#endif
