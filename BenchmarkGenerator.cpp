/*
Copyright (c) 2008-2010 Daniel Marbach & Thomas Schaffter

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper(s) listed
on http://gnw.sourceforge.net.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

$Id: BenchmarkGenerator.cpp 29 2012-01-04 17:06:55Z jbao $
*/

#include <ostream>
#include <sstream>
#include <string>
#include <omp.h>
#include "GnwSettings.h"
#include "BenchmarkGenerator.h"
#include "logging/logging.h"
using namespace ::logging;
	
// ============================================================================
// PUBLIC METHODS

/**
 * Constructor 
 */
BenchmarkGenerator::BenchmarkGenerator(const GeneNetwork& grn) {
	
	grn_ = grn;
	//steadyStateExperiments_ = new ArrayList<SteadyStateExperiment>();
	//steadyStateExperimentsODE_ = new ArrayList<SteadyStateExperiment>();
	//timeSeriesExperiments_ = new ArrayList<TimeSeriesExperiment>();
	//timeSeriesExperimentsODE_ = new ArrayList<TimeSeriesExperiment>();
	//wildTypeODE_ = NULL;
	//timeSeriesPerturbations_ = null;
}

BenchmarkGenerator::~BenchmarkGenerator() {

}

// ----------------------------------------------------------------------------

/** 
 * Run all experiments, save the gold standards and the datasets.
 * @throws CancelException, Exception 
 */
void BenchmarkGenerator::generateGoldStandard() {
	
	//GraphUtilities util = new GraphUtilities(grn_);
	//util.anonymizeGenes();
	
	::logging::log::emit<Info>() << "Starting benchmark generation ..." <<
		::logging::log::endl;
	
	//checkForInterruption();

	// save signed network
	std::string filename = GnwSettings::Instance()->getOutputDirectory() + grn_.getId() + "_goldstandard_signed.tsv";
	grn_.writeTSV(filename.c_str());
	
	//checkForInterruption();
	
	// save the complete network in smbl2
	filename = GnwSettings::Instance()->getOutputDirectory() + grn_.getId() + ".xml";
	grn_.writeSBML(filename.c_str());
	
	//checkForInterruption();
	
	// create and run the experiments
	// loadInitialConditions("tmp/InSilicoSize10-Yeast3-initial-conditions.tsv");
	runAll();
	
	//checkForInterruption();
	
	GnwSettings *set = GnwSettings::Instance();
	bool toAddExperimentalNoise = set->getAddNormalNoise() || set->getAddLognormalNoise() || set->getAddMicroarrayNoise(); 
	
	// print the data
	std::string postfix = "";
	if (toAddExperimentalNoise)
		postfix = "_noexpnoise";
	
	//checkForInterruption();
	
	for (unsigned int i=0; i<steadyStateExperiments_.size(); i++) {
		SteadyStateExperiment *exp = steadyStateExperiments_.at(i);
		exp->printAll(postfix);
		std::string label = exp->getLabel();
		//if (label == "multifactorial" || label == "dream4_timeseries" || label == "dualknockouts")
			//exp.getPerturbation().printPerturbations(exp.getLabel());
	}
	
	//checkForInterruption();
	
	//for (int i=0; i<steadyStateExperimentsODE_.size(); i++)
	//	steadyStateExperimentsODE_.get(i).printAll("_nonoise");

	//checkForInterruption();
	
    int pn = GnwSettings::Instance()->getPerturbationNumber();
	for (unsigned int i=0; i<timeSeriesExperiments_.size(); i++) {
		TimeSeriesExperiment exp = timeSeriesExperiments_.at(i);
        std::stringstream ss;
        ss << pn;
		exp.printAll(postfix + "_perturbation-" + ss.str());
		std::string label = exp.getLabel();
		//if (label == "dream4_timeseries" || label == "multifactorial_timeseries" || label == "dualknockout_timeseries")
			//exp.getPerturbation().printPerturbations(exp.getLabel());
	}
	
	//checkForInterruption();
	
	//for (int i=0; i<timeSeriesExperimentsODE_.size(); i++)
	//	timeSeriesExperimentsODE_.get(i).printAll("_nonoise");

	//checkForInterruption();
	
	// add noise, normalize, and print
	if (toAddExperimentalNoise) {
		
		addExperimentalNoise();
		if (set->getNormalizeAfterAddingNoise())
			normalize();
		
		for (unsigned int i=0; i<steadyStateExperiments_.size(); i++)
			steadyStateExperiments_.at(i)->printAll("");
		for (unsigned int i=0; i<timeSeriesExperiments_.size(); i++)
			timeSeriesExperiments_.at(i).printAll("");
	}
	
}

// ============================================================================
// PRIVATE METHODS

/** 
 * Run all experiments. If the timeSeriesExperiments_ have
 * already been created (e.g. with loadInitialConditions()), they are run as is.
 * Otherwise, they are created with random initial conditions according to the
 * settings in GnwSettings.
 * @throws CancelException, Exception 
 */
void BenchmarkGenerator::runAll() {

	//GnwSettings set = GnwSettings.getInstance();
	
	// The user should have selected either ODEs, SDEs, or both
	//if (!set.getSimulateODE() && !set.getSimulateSDE())
	//	throw new IllegalArgumentException("At least one of simulateODE_ and simulateSDE_ must be selected in GnwSettings");
	
	runSteadyStateExperiments();
	//checkForInterruption();
	runTimeSeriesExperiments();
}

// ----------------------------------------------------------------------------

/**
 * Create and run all steady-state experiments
 * @throws CancelException, Exception 
 */
void BenchmarkGenerator::runSteadyStateExperiments() {
	
	//GnwSettings set = GnwSettings.getInstance();
			
	// First, ODEs are always simulated (even if it's not set in the GnwSettings), because
	// we use the ODE wild-type as intial condition for the SDE wild-type, and the time-to-
	// convergence of the ODE as limit for the SDEs
	
	//checkForInterruption();
	
	// the wild-type
	std::string label = "wildtype";
	SteadyStateExperiment *wt = new SteadyStateExperiment(label);
	wt->setGrn(grn_);
	wt->run();
	Mat_DP perturb = wt->getSsPerturbation();
	constructInitialConditionFromWildType(perturb, wildTypeODE_);
	steadyStateExperiments_.push_back(wt);
	
    // Jacobian
    Vec_DP ss = wt->getSteadyStates();
    Vec_DP fvec(ss.size());
    grn_.computeDxydt(ss, fvec);
    Mat_DP df(ss.size(), ss.size());
    //wt->fdjac(ss, fvec, df, &GeneNetwork::dynamics);
    wt->jac(ss, df);
	std::string absPath = GnwSettings::Instance()->getOutputDirectory();
    wt->printJacobian(absPath + grn_.getId() + "_jacobian.tsv", df);

    //delete wt;
	/*		
	checkForInterruption();
	
	// knockouts
	if (set.generateSsKnockouts()) {
		PerturbationSingleGene knockouts = new PerturbationSingleGene(grn_);
		knockouts.singleGenePerturbations(0);
		createAndRunSsExperiment(Solver.type.ODE, knockouts, "knockouts");
	}
	
	checkForInterruption();
	
	// knockdowns
	if (set.generateSsKnockdowns()) {
		PerturbationSingleGene knockdowns = new PerturbationSingleGene(grn_);
		knockdowns.singleGenePerturbations(0.5);
		createAndRunSsExperiment(Solver.type.ODE, knockdowns, "knockdowns");
	}
	
	checkForInterruption();
	*/
	// multifactorial weak (DREAM4)
	if (GnwSettings::Instance()->generateSsMultifactorial()) {
		PerturbationMultifactorial *multifact = new PerturbationMultifactorial(&grn_);

		std::string label = "multifactorial";
		//if (!GnwSettings::Instance()->getLoadPerturbations())
			multifact->multifactorialAllGenesWeak(grn_.getSize());
		//else
			//multifact.loadPerturbations(label);
					
		//createAndRunSsExperiment(Solver.type.ODE, multifact, label);
		delete multifact;
	}
	
}

// ----------------------------------------------------------------------------

/**
 * Run all time-series experiments
 * @throws CancelException 
 */
void BenchmarkGenerator::runTimeSeriesExperiments() throw(std::exception) {
	
	//if (steadyStateExperiments_ == NULL)
	//	throw new RuntimeException("The wild-type must be simulated to run time-series experiments");
	
	// wild-type will be used as initial condition below (note, this is the SDE wild-type if SDEs are used)
	SteadyStateExperiment *wildType = steadyStateExperiments_.at(0);
	//if (wildType.getLabel() != "wildtype")
	//	throw new RuntimeException("The first steady-state experiment must be the wild-type");
	Vec_DP xy0;
	Mat_DP perturb = wildType->getSsPerturbation();
	constructInitialConditionFromWildType(perturb, xy0);
	/*
	GnwSettings set = GnwSettings.getInstance();
	Solver.type simulationType = Solver.type.ODE; 
	if (set.getSimulateSDE())
		simulationType = Solver.type.SDE;		
	
	checkForInterruption();
	
	// knockouts
	if (set.generateTsKnockouts()) {
		PerturbationSingleGene knockouts = new PerturbationSingleGene(grn_);
		knockouts.singleGenePerturbations(0);
		createAndRunTsExperiment(simulationType, knockouts, false, "knockout_timeseries", xy0);
	}
	
	checkForInterruption();
	
	// knockdowns
	if (set.generateTsKnockdowns()) {
		PerturbationSingleGene knockdowns = new PerturbationSingleGene(grn_);
		knockdowns.singleGenePerturbations(0.5);
		createAndRunTsExperiment(simulationType, knockdowns, false, "knockdown_timeseries", xy0);
	}
	
	checkForInterruption();
	*/
	// multifactorial weak
	if (GnwSettings::Instance()->generateTsMultifactorial()) {
		Perturbation *multifact = new PerturbationMultifactorial();
		std::string label = "multifactorial_timeseries";
		
		//if (!GnwSettings::Instance()->getLoadPerturbations()) {
			// use the same perturbations as for the steady-state experiments, if they were simulated
			for (unsigned int i=0; i<steadyStateExperiments_.size(); i++) {
				
				//checkForInterruption();
				
				if (steadyStateExperiments_.at(i)->getLabel() == "multifactorial") {
					multifact = steadyStateExperiments_.at(i)->getPerturbation();
					break;
				}
			}
			if (multifact->getNumGenes() == 0) {
				multifact = new PerturbationMultifactorial(&grn_);
				//multifact->multifactorialHub(GnwSettings::Instance()->getNumTimeSeries());
				multifact->multifactorialAllGenesWeak(GnwSettings::Instance()->getNumTimeSeries());
				//multifact->perturbSingleGene(GnwSettings::Instance()->getNumTimeSeries(), "crp");
    
                int pn = GnwSettings::Instance()->getPerturbationNumber();
                std::stringstream ss;
                ss << pn;
                multifact->printPerturbations("perturbation-" + ss.str() + "_" + 
                        label);
			}
		//} else {
		//	multifact = new PerturbationMultifactorial(grn_);
			//multifact.loadPerturbations("multifactorial");
		//}
		//multifact->printPerturbations("");
		createAndRunTsExperiment(multifact, false, label, xy0);
		delete multifact;
	}
	
	//checkForInterruption();
	/*
	// multifactorial strong
	if (GnwSettings::Instance()->generateTsDREAM4TimeSeries()) {
		PerturbationMultifactorial *multifact = NULL;
		std::string label = "dream4_timeseries";
		
		if (!GnwSettings::Instance()->getLoadPerturbations()) {
			// use the same perturbations as for the steady-state experiments, if they were simulated
			for (int i=0; i<steadyStateExperiments_.size(); i++) {
				
				//checkForInterruption();
				
				if (steadyStateExperiments_.at(i).getLabel() == "dream4_timeseries") {
					multifact = &(steadyStateExperiments_.at(i).getPerturbation());
					break;
				}
			}
			if (multifact == NULL) {
				multifact = new PerturbationMultifactorial(grn_);
				multifact->multifactorialStrong(GnwSettings::Instance()->getNumTimeSeries());
			}
		} else {
			multifact = new PerturbationMultifactorial(grn_);
			//multifact.loadPerturbations(label);
		}
		createAndRunTsExperiment(*multifact, true, label, xy0);
	}
	*/
}

// ----------------------------------------------------------------------------

/**
 * The given experiment should be the wild-type, returns the concatenated mRNA
 * and protein concentrations
 */
void BenchmarkGenerator::constructInitialConditionFromWildType(Mat_I_DP& perturb, Vec_O_DP& xy0) {
	
	//Mat_DP perturb = wildType.getSsPerturbation();
	Vec_DP x(perturb.ncols());
	for (int i = 0; i < x.size(); ++i) {
		x[i] = perturb[0][i];
	}

	if (GnwSettings::Instance()->getModelTranslation()) {
		
		//Mat_DP perturb = wildType.getSsPerturbationProteins();
		Vec_DP y(perturb.ncols());
		for (int i = 0; i < y.size(); ++i) {
			y[i] = perturb[0][i];
		}
		
		xy0 = Vec_DP(2*grn_.getSize());
		
		for (int i=0; i<x.size(); i++)
			xy0[i] = x[i];
		
		for (int i=0; i<y.size(); i++)
			xy0[x.size()+i] = y[i];
		
		//return xy0;
		
	} else {
		//xy0 = Vec_DP(grn_.getSize());
		xy0 = x;
	}
}

// ----------------------------------------------------------------------------

/** Create and run a time-series experiment, add it to timeSeriesExperiments_ */
void BenchmarkGenerator::createAndRunTsExperiment(Perturbation *perturbation, bool restoreWildTypeAtHalftime, std::string label, Vec_DP xy0) {
    int nts = GnwSettings::Instance()->getNumTimeSeries();
    int pn = GnwSettings::Instance()->getPerturbationNumber();
    //TimeSeriesExperiment ts;
#pragma omp parallel for num_threads(nts) 
	for (int i=0; i<nts; i++) {
	    TimeSeriesExperiment ts = TimeSeriesExperiment(perturbation, restoreWildTypeAtHalftime, label);
        GeneNetwork grn = GeneNetwork(grn_);
	    ts.setGrn(grn);
	    
        std::stringstream ss, ssrank;
		//ss.str("");
		ss << i + 1;
        int rank = omp_get_thread_num();
        ssrank << rank;
		::logging::log::emit<Info>() << "Rank " << ssrank.str().c_str() << " simulating time-series number " << ss.str().c_str() << " ..." << ::logging::log::endl;
	    ts.run(xy0);
        ts.integrate(i);
#pragma omp critical
	    timeSeriesExperiments_.push_back(ts);
        //delete ts;
    }
	//delete ts;
}

// ----------------------------------------------------------------------------

/**
 * Add log-normal noise to the data
 */
void BenchmarkGenerator::addExperimentalNoise() {

	for (unsigned int i=0; i<steadyStateExperiments_.size(); i++)
		steadyStateExperiments_.at(i)->addNoise();

	for (unsigned int i=0; i<timeSeriesExperiments_.size(); i++)
		timeSeriesExperiments_.at(i).addNoise();
}

// ----------------------------------------------------------------------------

/**
 * Normalize the data (divide by the maximum)
 */
void BenchmarkGenerator::normalize() {

	double max = 0.;

	// find the maximum concentration value off all experiments
	for (unsigned int i=0; i<steadyStateExperiments_.size(); i++) {
		double max_i = steadyStateExperiments_.at(i)->getMaximumConcentration();
		if (max_i > max)
			max = max_i;
	}

	for (unsigned int i=0; i<timeSeriesExperiments_.size(); i++) {
		double max_i = timeSeriesExperiments_.at(i).getMaximumConcentration();
		if (max_i > max)
			max = max_i;
	}

	std::stringstream ss;
	ss << max;
	::logging::log::emit<Info>() << "Normalizing with respect to max = " << ss.str().c_str() << ::logging::log::endl;
	
	// save the coefficient
	std::string filename = GnwSettings::Instance()->getOutputDirectory() + grn_.getId() + "_normalization_constant.tsv";
	::logging::log::emit<Info>() << "Writing file " << filename.c_str() << ::logging::log::endl;
	
	try {
		ofstream data_file(filename.c_str()); 
		data_file << max << std::endl;
		data_file.close();
	} catch (Exception& e) {
		::logging::log::emit<Info>() << "Error writing file, exception " << e.what().c_str() << ::logging::log::endl;
		throw Exception();
	}

	// normalize according to this max
	for (unsigned int i=0; i<steadyStateExperiments_.size(); i++)
		steadyStateExperiments_.at(i)->normalize(max);
	//for (int i=0; i<steadyStateExperimentsODE_.size(); i++)
	//	steadyStateExperimentsODE_.at(i).normalize(max);

	if (timeSeriesExperiments_.size() != 0) {
		for (unsigned int i=0; i<timeSeriesExperiments_.size(); i++)
			timeSeriesExperiments_.at(i).normalize(max);
	}
	
}
