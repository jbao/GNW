// $Id: SteadyStateExperiment.h 21 2011-08-08 12:08:06Z jbao $

#ifndef STEADYSTATEEXPERIMENT_H
#define STEADYSTATEEXPERIMENT_H

#include "nr/nr.h"
#include "GnwSettings.h"
#include "Exception.h"
#include "Experiment.h"

/** 
 * Implements experiments where the stead-state of the network is measured.
 * @author Daniel Marbach (firstname.name@gmail.com)
 * @author Thomas Schaffter (firstname.name@gmail.com)
 * 
 */
 class SteadyStateExperiment : public Experiment {
 
 public:
 	SteadyStateExperiment(const std::string& label);
 	SteadyStateExperiment(PerturbationMultifactorial *perturbation, std::string label);
 	~SteadyStateExperiment() {}
 
 	Mat_DP getSsPerturbation() { return ssPerturbation_; }
	Mat_DP getSsPerturbationProteins() { return ssPerturbationProteins_; }
	std::vector<double>& getTimeToConvergenceODE() { return timeToConvergenceODE_; }
	void setTimeToConvergenceODE(const std::vector<double>& t) { timeToConvergenceODE_ = t; }
	
	void run();
	void run(Vec_DP& xy0);
	void printAll(std::string postfix);
	void printMRNA(std::string postfix);
	void printProteins(std::string postfix);
	void computeSteadyStates();
	void constructInitialCondition(Vec_DP& xy0);
	void printSteadyStates(std::string filename, const Vec_DP& wt);
	void printSteadyStates(std::string filename, const Mat_DP& data);
	void printJacobian(std::string filename, const Mat_DP& data);
	
	bool converged(Vec_DP& previousState, Vec_DP& state);
	
	double getMaximumConcentration();
	void addNoise();
	void normalize(double max);
 
    void fdjac(Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df, 
            void (*dynamics)(Vec_I_DP&,Vec_O_DP&,GeneNetwork&));
    void jac(Vec_I_DP&, Mat_O_DP&);

    Vec_DP& getSteadyStates() { return steadyStates_; }
 
private:
 	/** Perturbed steady-states */
	Mat_DP ssPerturbation_;
	/** Perturbed steady-states for the proteins */
	Mat_DP ssPerturbationProteins_;
    /** Save the steady state concentrations */
    Vec_DP steadyStates_;

	/** Time of the current steady-state computation */
	double t_;
	/**
	 * For SteadyStateExperimentODE: return the steady-states as soon as convergence is reached.
	 * If there is no convergence until time maxt_, the values at this point are returned and a
	 * warning message is displayed.
	 */
	double maxtODE_;
	/**
	 * For stochastic simulations (SDEs), If maxtSDE < 0, we return the state at time 1.5*timeToConvergenceODE_,
	 * where timeToConvergenceODE_ should be set to the time of convergence for the deterministic simulation 
	 * of the same experiment. If maxtSteadyStateSDE > 0, we return the state at that time.
	 */
	//private double maxtSDE_;
	/** 
	 * For ODEs: save the time until convergence for each perturbation
	 * For SDEs: return the state at these times.
	 */
	std::vector<double> timeToConvergenceODE_;
	
	void computeSteadyState(double maxt);
 
 };
 
 #endif
