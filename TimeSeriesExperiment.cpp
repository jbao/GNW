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

$Id: TimeSeriesExperiment.cpp 29 2012-01-04 17:06:55Z jbao $
*/

#include <assert.h>
#include <fstream>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <omp.h>
#include "GnwSettings.h"
#include "Exception.h"
#include "TimeSeriesExperiment.h"
#include "logging/logging.h"
using namespace ::logging;
	
//extern "C"
inline int gene_network_dynamics(double, const double xy[], double dxydt[], void *param) { 
	
	TimeSeriesExperiment *exp = (TimeSeriesExperiment*) param;
	GeneNetwork* grn = &(exp->getGrn());
  	bool modelTranslation = GnwSettings::Instance()->getModelTranslation();
	//double[] dxydt = new double[xy.length];
	int size = grn->getSize();
	
    Vec_DP x(size), y(size);
	for (int i=0; i<size; i++) {
		if (xy[i] >= 0)
            x[i] = xy[i];
        else
            x[i] = 0;
		//std::stringstream s1,s2;
		//s1 << i;
		//s2 << x_[i];
		//::logging::log::emit<Debug>() << "x" << s1.str().c_str() << " = " << s2.str().c_str() << ::logging::log::endl;
	}
	grn->setX(x);

	if (modelTranslation)
		for (int i=0; i<size; i++)
			y[i] = xy[size+i];
	else
		y = x;
    grn->setY(y);
	
	// dxydt temporarily used to store the production rates of mRNA
	Vec_DP dxydt_v(size); 
	grn->computeMRnaProductionRates(dxydt_v);
	
	for (int i=0; i < size; i++)
		dxydt_v[i] = dxydt_v[i] - grn->getNodes().at(i).computeMRnaDegradationRate(grn->getX()[i]);
		
	for (int i=0; i < size; i++)
		dxydt[i] = dxydt_v[i];
	
	/*
	if (modelTranslation)
		for (int i=0; i<size; i++)
			dxydt[size+i] = nodes_.at(i).getMaxTranslation()*x_[i] - nodes_.at(i).computeProteinDegradationRate(y_[i]);
  	*/
  	return GSL_SUCCESS;
}
	
// ----------------------------------------------------------------------------

inline int gene_network_dynamics_constant_input(double, const double xy[], double dxydt[], void *param) { 
	
	TimeSeriesExperiment *exp = (TimeSeriesExperiment*) param;
    Perturbation* perturb = exp->getPerturbation();
	GeneNetwork* grn = &(exp->getGrn());
  	bool modelTranslation = GnwSettings::Instance()->getModelTranslation();
	//double[] dxydt = new double[xy.length];
	int size = grn->getSize();
	
    Vec_DP x(size), y(size);
	for (int i=0; i<size; i++) {
        if (xy[i] >= 0)
		    x[i] = xy[i];
        else
            x[i] = 0;
		//std::stringstream s1,s2;
		//s1 << i;
		//s2 << x_[i];
		//::logging::log::emit<Debug>() << "x" << s1.str().c_str() << " = " << s2.str().c_str() << ::logging::log::endl;
	}
	grn->setX(x);

	if (modelTranslation)
		for (int i=0; i<size; i++)
			y[i] = xy[size+i];
	else
		y = x;
	grn->setY(y);

	// dxydt temporarily used to store the production rates of mRNA
	Vec_DP dxydt_v(size); 
	grn->computeMRnaProductionRates(dxydt_v);
	
	for (int i=0; i < size; i++)
		dxydt_v[i] = dxydt_v[i] - grn->getNodes().at(i).computeMRnaDegradationRate(grn->getX()[i]) + perturb->getPerturbations()[0][i];
		
	for (int i=0; i < size; i++)
		dxydt[i] = dxydt_v[i];
	
	/*
	if (modelTranslation)
		for (int i=0; i<size; i++)
			dxydt[size+i] = nodes_.at(i).getMaxTranslation()*x_[i] - nodes_.at(i).computeProteinDegradationRate(y_[i]);
  	*/
  	return GSL_SUCCESS;
}
	
// ============================================================================
// PUBLIC METHODS

/**
 * Constructor
 */
TimeSeriesExperiment::TimeSeriesExperiment(Perturbation *perturbation, bool restoreWildTypeAtHalftime, std::string label) : Experiment(perturbation, label) {
	
	//timeSeries_ = null;
	//timeSeriesProteins_ = null;
	//xy0_ = null;
	restoreWildTypeAtHalftime_ = restoreWildTypeAtHalftime;
	setMaxtAndNumTimePoints();
}
	
TimeSeriesExperiment::~TimeSeriesExperiment() {

}
	
// ----------------------------------------------------------------------------

TimeSeriesExperiment::TimeSeriesExperiment(const TimeSeriesExperiment& e) : Experiment(e) {
    //grn_ = e.grn_;
    timeSeries_ = e.timeSeries_;
    inputs_ = e.inputs_;
    outputs_ = e.outputs_;
    numTimePoints_ = e.numTimePoints_;
}

// ----------------------------------------------------------------------------

TimeSeriesExperiment& TimeSeriesExperiment::operator=(const TimeSeriesExperiment& rhs) {
    Experiment::operator=(rhs);
    //grn_ = rhs.grn_;
    timeSeries_ = rhs.timeSeries_;
    inputs_ = rhs.inputs_;
    outputs_ = rhs.outputs_;
    numTimePoints_ = rhs.numTimePoints_;
}

// ----------------------------------------------------------------------------

/**
 * Run all experiments
 */
void TimeSeriesExperiment::run(Vec_DP xy0) {
	
	xy0_ = xy0;
	
	std::string simulationType = "ODEs";
	//if (solverType_ == Solver.type.SDE)
	//	simulationType = "SDEs";
	::logging::log::emit<Info>() << "Simulating time-series " << label_.c_str() << " using " << simulationType.c_str() << " ..." << ::logging::log::endl;

	//bool simulateLoadedExperiments = (&timeSeries_ != NULL);
	//if (simulateLoadedExperiments)
	//	throw new RuntimeException("NEEDS TO BE FIXED, NOT FUNCTIONAL");
	
	//if (!simulateLoadedExperiments) {
		//timeSeries_ = new ArrayList<DoubleMatrix2D>();
	//	if (GnwSettings::Instance()->getModelTranslation())
	//		timeSeriesProteins_ = new ArrayList<DoubleMatrix2D>();
	//}
	
	// create and run the time series experiments
    //int i;
//#pragma omp parallel for num_threads(numExperiments_) 
	//for (i=0; i<numExperiments_; i++) {
	//    std::stringstream ss, ssrank;
		//ss.str("");
	//	ss << i + 1;
    //    int rank = omp_get_thread_num();
    //    ssrank << rank;
	//	::logging::log::emit<Info>() << "Rank " << ssrank.str().c_str() << " simulating time-series number " << ss.str().c_str() << " ..." << ::logging::log::endl;
//#pragma omp critical
    //    integrate(k);
	//}
	//log.log(Level.INFO, "");
	
}
		
// ----------------------------------------------------------------------------

/**
 * Run the numerical integration of the k'th time-series and add the results to timeSeries_ and timeSeriesProteins_.
 * The wild-type is restored after the experiments.
 */
void TimeSeriesExperiment::integrate(int k) {

	if (GnwSettings::Instance()->getDt()*(numTimePoints_-1) != maxt_) {
		::logging::log::emit<Error>() << "dt * (numTimePoints-1) != maxt" << ::logging::log::endl;
		exit(1);
	}
	
	// allocate space
	Mat_DP ts(numTimePoints_, numGenes_);
	Mat_DP input(numTimePoints_, numGenes_);
	Mat_DP output(numTimePoints_, numGenes_);
	Mat_DP tsProteins;
	//DoubleMatrix2D tsProteins = null;
	if (GnwSettings::Instance()->getModelTranslation())
		Mat_DP tsProteins(numTimePoints_, numGenes_);

	if (xy0_.size() == 0) 
		throw Exception("TimeSeriesExperiment:integrate(): No initial condition set!");
	
	if (maxt_ <= 0) {
		::logging::log::emit<Info>() << "Duration (maxt) must be greater than 0! maxt_ = " << maxt_ << ::logging::log::endl;
		throw Exception();
	}

	//Solver solver = new Solver(solverType_, grn_, xy0_.toArray());
	double t = 0;
	
	/*
	// for SDEs, simulate the wild-type for a short time to get a new independent sample
	if (solverType_ == Solver.type.SDE) {
		double tlim = maxt_/10.0;
		do {
			try {
				t += solver.step();
			} catch (Exception e) {
				log.log(Level.INFO, "TimeSeriesExperiment.integrate(): Exception in phase 0, t = " + t + ":" + e.getMessage());
				throw new RuntimeException();
			}
		} while (t < tlim);

		// set this sample as the new initial condition
		xy0_.assign(solver.getState()); 
	}
	*/
	
	// Set first line of the time series dataset (at t=0)
	for (int i=0; i<numGenes_; i++) {
		xy0_[i] += perturbation_->getPerturbations()[k][i];
        //xy0_[i] += 1;
        if (xy0_[i] < 0)
            xy0_[i] = 0;
		ts[0][i] = xy0_[i];
    }
	if (GnwSettings::Instance()->getModelTranslation())
		for (int i=0; i<numGenes_; i++)
			tsProteins[0][i] = xy0_[numGenes_+i];
	
	// apply perturbation
	//perturbation_->applyPerturbation(k);
	//for (int i=0; i<numGenes_; i++)
	//	grn_.getNodes().at(i).perturbBasalActivation( 
    //            perturbation_->getPerturbations()[k][i] );
	
    t = 0; // reset time, the time-series only really starts here
	double dt = GnwSettings::Instance()->getDt();
	
	if (dt <= 0 || dt > maxt_) {
		::logging::log::emit<Error>() << "Interval between two measuread points must be >0 and <maxt." << ::logging::log::endl;
		exit(1);
	}
	
	// TODO: verify why GnwSettings->dt_ is used => most likely to be removed
	// replace that by being sure that dt is a multiple of the internal integration step-size
	// as dt is calculated from maxt_ and numTimePoints there must be no problem (see reported bug)
	double frac = maxt_/dt;
	if (frac - static_cast<int>(frac) != 0) {
		::logging::log::emit<Error>() << "Duration (maxt) must be a multiple of numTimePoints-1." << ::logging::log::endl;
		exit(1);
	}
	
	//double tlim = maxt_/2.0 - 1e-12;
	bool wildTypeRestored = false;
	int pt = 1;
	
	/*
	do {
		double t1 = t;
		try {
			// For ODEs: this steps the time by dt_, but using an adaptive internal step size
			// to guarantee the specified tolerance (getRate() may be called several times for one step)
			// For SDEs: this steps the time by dt_, the solver integrates with a smaller, fixed step size
			// defined in SDESettings by dt_*multiplier_ (SDESettings.dt_ != TimeSeriesExperiment.dt_)
			t += solver.step();
		} catch (Exception e) {
			log.log(Level.INFO, "TimeSeriesExperiment.integrate(): Exception at t = " + t + ":" + e.getMessage());
			throw new RuntimeException();
		}
		
		if (t != t1 + dt)
			throw new RuntimeException("Solver failed to step time by dt, expected t = " + (t1+dt) + ", obtained t = " + t);
		
		if (restoreWildTypeAtHalftime_ && t >= tlim && !wildTypeRestored) {
			perturbation_.restoreWildType();
			wildTypeRestored = true;
		}
		
		// Save the state of the result
		double[] xy = solver.getState();
		for (int g=0; g<numGenes_; g++)
			ts.set(pt, g, xy[g]);

		if (modelTranslation_)
			for (int g=0; g<numGenes_; g++)
				tsProteins.set(pt, g, xy[numGenes_+g]);
		
		pt++;
	} while (t < maxt_);
	*/
	/*
	// Euler
	for (; t < maxt_; t += dt) {
	
		if (restoreWildTypeAtHalftime_ && t >= tlim && !wildTypeRestored) {
			perturbation_.restoreWildType();
			wildTypeRestored = true;
		}
		
		// Save the state of the result
		Vec_DP state(numGenes_), rate(numGenes_);
		for (int i = 0; i < numGenes_; ++i)
			state[i] = ts[pt-1][i];
		grn_.computeDxydt(state, rate);
		//double[] xy = solver.getState();
		for (int g=0; g<numGenes_; g++)
			ts[pt][g] = state[g] + rate[g]*dt;

		//if (GnwSettings::Instance()->getModelTranslation())
		//	for (int g=0; g<numGenes_; g++)
		//		tsProteins[pt][g] = xy[numGenes_+g];
		
		pt++;
	}
	
	// rkf45
	const int NVAR=numGenes_, NSTEP=frac;
    int j;
    DP x1=t, x2=maxt_;
    Vec_DP vstart=xy0_;

    // Note: The arrays xx and y must have indices up to NSTEP
    xx_p = new Vec_DP(NSTEP+1);
    y_p = new Mat_DP(NVAR,NSTEP+1);
    Vec_DP &xx = *xx_p;
    Mat_DP &y = *y_p;
    void (GeneNetwork::*funcPtr)( const DP, Vec_I_DP&, Vec_O_DP& ) = &GeneNetwork::dynamics; 
    Experiment::rkdumb(vstart,x1,x2,funcPtr);
    
    for (j=1; j <= NSTEP; j++) {
        for (int g=0; g<numGenes_; g++)
			ts[j][g] = y[g][j];
    }
    delete y_p;
    delete xx_p;
	*/
	// gsl
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
     
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, numGenes_);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (GnwSettings::Instance()->getAbsolutePrecision(), GnwSettings::Instance()->getRelativePrecision());
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (numGenes_);
     
    gsl_odeiv_system sys;
	if (!GnwSettings::Instance()->generateTsConstantInput()) 
        sys = {gene_network_dynamics, NULL, numGenes_, this};
    else
        sys = {gene_network_dynamics_constant_input, NULL, numGenes_, this};
        
    double y[numGenes_];
    Vec_DP y_vec(numGenes_);
    for (int i=0; i<numGenes_; i++) {
		y[i] = xy0_[i];
		y_vec[i] = xy0_[i];
        std::stringstream s1, s2;
        s1 << i;
        s2 << y[i];
		//::logging::log::emit<Debug>() << "y" << s1.str().c_str() << " = " << s2.str().c_str() << ::logging::log::endl;
	}
	/*
    while (t < maxt_) {
        int status = gsl_odeiv_evolve_apply (e, c, s,
                                                &sys, 
                                                &t, maxt_,
                                                &dt, y);
     
        if (status != GSL_SUCCESS)
            break;
            
        std::stringstream ss;
        ss << t;
        ::logging::log::emit<Debug>() << "t = " << ss.str().c_str() << ::logging::log::endl;
        
        for (int i=0; i<numGenes_; i++)
        	ts[pt][i] = y[i];
       	pt++;
     
    }
    */
    // non-adaptive steps
    double y_err[numGenes_];
    double dydt_in[numGenes_], dydt_out[numGenes_];
     
    /* initialise dydt_in from system parameters */
    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    // initialize inputs
    Vec_DP current_input(numGenes_);
    grn_.computeInputs(y_vec, current_input);
    for (int i=0; i<numGenes_; i++) {
		input[0][i] = current_input[i];
        int numOutputs = grn_.getNodes()[i].getNumOutputs();
        output[0][i] = y_vec[i] * numOutputs;
    }
     
    while (t < maxt_) {
        int status = gsl_odeiv_step_apply (s, t, dt, 
                                              y, y_err, 
                                              dydt_in, 
                                              dydt_out, 
                                              &sys);
     
        if (status != GSL_SUCCESS)
            break;
     
     	std::stringstream ss;
        ss << t;
        //::logging::log::emit<Debug>() << "t = " << ss.str().c_str() << ::logging::log::endl;
     
     	// save results
     	for (int i=0; i<numGenes_; i++) {
            if (y[i] < 0)
                y[i] = 0;
			y_vec[i] = y[i];
        }
     	grn_.computeInputs(y_vec, current_input);
     	for (int i=0; i<numGenes_; i++) {
     		ts[pt][i] = y[i];
        	dydt_in[i] = dydt_out[i];
        	input[pt][i] = current_input[i];
            int numOutputs = grn_.getNodes()[i].getNumOutputs();
            output[pt][i] = y_vec[i] * numOutputs;
        }
     
        t += dt;
        pt++;
     
    }

     
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);


	//assert(t == maxt_); 
	//assert(pt == numTimePoints_);
	
	// make sure the wild-type is restored
	//if (!wildTypeRestored)
	//	perturbation_->restoreWildType();
	
	// add the new time-series data to the array lists
	timeSeries_.push_back(ts);
	inputs_.push_back(input);
	outputs_.push_back(output);
	if (GnwSettings::Instance()->getModelTranslation())
		timeSeriesProteins_.push_back(tsProteins);
	
	//getSDESolver()
	//if (solverType_ == Solver.type.SDE && solver.getXNegativeCounter() > 0)
	//	log.log(Level.INFO, "SDE: " + solver.getXNegativeCounter() + " times a concentration became negative due to noise and was set to 0");

}
	
// ----------------------------------------------------------------------------

/**
 * Run the numerical integration of the k'th time-series and add the results to timeSeries_ and timeSeriesProteins_.
 * The wild-type is restored after the experiments.
 */
//void TimeSeriesExperiment::integrateConstantInput(int k) {
//
//	// allocate space
//	//Mat_DP ts(numTimePoints_, numGenes_);
//	//Mat_DP input(numTimePoints_, numGenes_);
//	//Mat_DP output(numTimePoints_, numGenes_);
//	//Mat_DP tsProteins;
//	////DoubleMatrix2D tsProteins = null;
//	//if (GnwSettings::Instance()->getModelTranslation())
//	//	Mat_DP tsProteins(numTimePoints_, numGenes_);
//    std::vector<Vec_DP> ts_vec;
//
//	if (xy0_.size() == 0) 
//		throw Exception("TimeSeriesExperiment:integrateConstantInput(): No initial condition set!");
//	
//	double t = 0;
//	double dt = GnwSettings::Instance()->getDt();
//	
//	int pt = 1;
//	
//	// gsl
//	const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
//     
//    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, numGenes_);
//    gsl_odeiv_control * c = gsl_odeiv_control_y_new (GnwSettings::Instance()->getAbsolutePrecision(), GnwSettings::Instance()->getRelativePrecision());
//    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (numGenes_);
//     
//    gsl_odeiv_system sys = {GeneNetwork::gene_network_dynamics, NULL, numGenes_, &grn_};
//        
//    double y[numGenes_];
//    Vec_DP y_vec(numGenes_);
//    for (int i=0; i<numGenes_; i++) {
//		y[i] = xy0_[i];
//		y_vec[i] = xy0_[i];
//        std::stringstream s1, s2;
//        s1 << i;
//        s2 << y[i];
//		//::logging::log::emit<Debug>() << "y" << s1.str().c_str() << " = " << s2.str().c_str() << ::logging::log::endl;
//	}
//    // non-adaptive steps
//    double y_err[numGenes_];
//    double dydt_in[numGenes_], dydt_out[numGenes_];
//     
//    /* initialise dydt_in from system parameters */
//    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
//    // initialize inputs
//    Vec_DP current_input(numGenes_);
//    grn_.computeInputs(y_vec, current_input);
//    //for (int i=0; i<numGenes_; i++) {
//	//	input[0][i] = current_input[i];
//    //    int numOutputs = grn_.getNodes()[i].getNumOutputs();
//    //    output[0][i] = y_vec[i] * numOutputs;
//    //}
//     
//    Vec_DP previous_y;
//    ts_vec.push_back(y_vec);
//    do {
//        
//        for (int i=0; i<numGenes_; i++) {
//            y[i] += perturbation_->getPerturbations()[k][i];
//            y_vec[i] = y[i];
//        }
//	
//        previous_y = y_vec;
//        
//        int status = gsl_odeiv_step_apply (s, t, dt, 
//                                              y, y_err, 
//                                              dydt_in, 
//                                              dydt_out, 
//                                              &sys);
//     
//        if (status != GSL_SUCCESS)
//            break;
//     
//     	std::stringstream ss;
//        ss << t;
//        //::logging::log::emit<Debug>() << "t = " << ss.str().c_str() << ::logging::log::endl;
//     
//     	// save results
//     	for (int i=0; i<numGenes_; i++) {
//            if (y[i] < 0)
//                y[i] = 0;
//			y_vec[i] = y[i];
//        }
//     	grn_.computeInputs(y_vec, current_input);
//     	for (int i=0; i<numGenes_; i++) {
//     		//ts[pt][i] = y[i];
//        	dydt_in[i] = dydt_out[i];
//        	//input[pt][i] = current_input[i];
//            //int numOutputs = grn_.getNodes()[i].getNumOutputs();
//            //output[pt][i] = y_vec[i] * numOutputs;
//        }
//     
//        t += dt;
//        pt++;
//     
//        ts_vec.push_back(y_vec);
//
//    } while(!converged(previous_y,y_vec));
//
//    gsl_odeiv_evolve_free (e);
//    gsl_odeiv_control_free (c);
//    gsl_odeiv_step_free (s);
//
//	// add the new time-series data to the array lists
//	timeSeries_.push_back(ts);
//	inputs_.push_back(input);
//	outputs_.push_back(output);
//	if (GnwSettings::Instance()->getModelTranslation())
//		timeSeriesProteins_.push_back(tsProteins);
//	
//}
	
// ----------------------------------------------------------------------------

/** 
 * Implements the method gsl_multiroot_test_delta() of GSL:
 * This function tests for the convergence of the sequence by comparing the last step dx
 * with the absolute error epsabs and relative error epsrel to the current position x.
 * The test returns true if the following condition is achieved:
 * 		|dx_i| < epsabs + epsrel |x_i|
 * for each component of x and returns false otherwise.
 */
//bool TimeSeriesExperiment::converged(Vec_DP& previousState, Vec_DP& state) {
//	 
//	double absolutePrecision = GnwSettings::Instance()->getAbsolutePrecision();
//	double relativePrecision = GnwSettings::Instance()->getRelativePrecision();
//	 
//	for (int i=0; i<previousState.size(); i++) {
//		
//		double dxy = abs(previousState[i] - state[i]); 
//		
//		if (dxy > absolutePrecision + relativePrecision*abs(state[i])) {
//			// remember point
//			for (int k=0; k<previousState.size(); k++)
//				previousState[k] = state[k];
//			
//			return false;
//		}
//	}
//	return true;
//}
	
// ----------------------------------------------------------------------------

/** 
 * Print all the trajectories to a single file, the initial conditions are printed
 * to a separate file. Protein trajectories are only printed if translation is
 * modelled. Append the given string to the filenames (e.g. "-nonoise"). 
 */
void TimeSeriesExperiment::printAll(std::string postfix) {
	
	if (timeSeries_.size() < 1)
		return;
	
	if (postfix != "_noexpnoise") {
		printTrajectories(postfix + "_" + label_, timeSeries_);    // print mRNA time courses
		printMatrix(postfix + "_inputs", inputs_);
		printMatrix(postfix + "_outputs", outputs_);
	}
	
	if (GnwSettings::Instance()->getModelTranslation())
		printTrajectories(postfix + "_proteins_" + label_, timeSeriesProteins_); // print protein time courses
}


// ----------------------------------------------------------------------------

/** 
 * Print all the trajectories to a single file, the initial conditions are printed
 * to a separate file. If the argument is set true, the protein instead of the
 * mRNA concentrations are printed. append the given string to the filenames (e.g. "-nonoise"). 
 */
void TimeSeriesExperiment::printTrajectories(std::string postfix, std::vector<Mat_DP> timeSeries) {
			
	try { 
		// Filename
		std::string filename = GnwSettings::Instance()->getOutputDirectory() + grn_.getId() + postfix + ".tsv";
		ofstream data_file(filename.c_str()); 
		
		// Header
		data_file << "\"Time\"\t" << grn_.getHeader(false);

		// For every time series...
		for (unsigned int i=0; i<timeSeries.size(); i++) {

			// The data
			Mat_DP data = timeSeries.at(i);
			double dt = GnwSettings::Instance()->getDt();

			data_file << std::endl;
			for (int tp=0; tp<numTimePoints_; tp++) {
				data_file << tp*dt;

				for (int g=0; g<numGenes_; g++) {
                    std::stringstream sstp, ssg;
                    sstp << tp;
                    ssg << g;
                    //::logging::log::emit<Debug>() << "tp = " << sstp.str().c_str() << " g = " << ssg.str().c_str() << ::logging::log::endl;
					data_file << "\t" << data[tp][g];
                }
				data_file << std::endl;
			}
		}

		data_file.close();
		::logging::log::emit<Info>() << "Writing file " << filename.c_str() << ::logging::log::endl;

	} catch (exception& e) {
		::logging::log::emit<Info>() << "TimeSeriesExperiment:printTrajectories(): " << e.what() << ::logging::log::endl;
		throw Exception();
	}
}

//	 ----------------------------------------------------------------------------

/** 
 * Print the input/output matrix over time to a single file.
 * @param postfix String appended to the filename
 * @param matrix Input/output matrix to be printed
 */
void TimeSeriesExperiment::printMatrix(std::string postfix, 
                                       std::vector<Mat_DP> matrix) {
			
	try { 
		// Filename
		std::string filename = GnwSettings::Instance()->getOutputDirectory() + grn_.getId() + postfix + ".tsv";
		ofstream data_file(filename.c_str()); 
		
		// Header
		data_file << "\"Time\"\t" << grn_.getHeader(false);

		// For every time series...
		for (unsigned int i=0; i<matrix.size(); i++) {

			// The data
			Mat_DP data = matrix.at(i);
			double dt = GnwSettings::Instance()->getDt();

			data_file << std::endl;
			for (int tp=0; tp<numTimePoints_; tp++) {
				data_file << tp*dt;

				for (int g=0; g<numGenes_; g++)
					data_file << "\t" << data[tp][g];
				data_file << std::endl;
			}
		}

		data_file.close();
		::logging::log::emit<Info>() << "Writing file " << filename.c_str() << ::logging::log::endl;

	} catch (exception& e) {
		::logging::log::emit<Info>() << "TimeSeriesExperiment:printMatrix(): " << e.what() << ::logging::log::endl;
		throw Exception();
	}
}

// ============================================================================
// PRIVATE FUNCTIONS

/** Set maxt_ and numTimePoints_ according to GnwSettings (checks that they are consistent with GnwSettings.dt_) */
void TimeSeriesExperiment::setMaxtAndNumTimePoints() {
	
	double dt = GnwSettings::Instance()->getDt();
	maxt_ = GnwSettings::Instance()->getMaxtTimeSeries();
	numTimePoints_ = static_cast<int>(maxt_/dt) + 1;

	if (dt*(numTimePoints_-1) != maxt_) {
		::logging::log::emit<Error>() << "maxt must be a multiple of dt" << ::logging::log::endl;
		exit(1);
	}
}
	
// ----------------------------------------------------------------------------

/**
 * Get the maximum concentration in the time series.
 * For now, only get the max concentration between mRNA levels. Later, perhaps
 * leave the choice to the user if he prefer, e.g., to normalise by the mRNA
 * OR¬†protein levels.
 */
double TimeSeriesExperiment::getMaximumConcentration() {

	double max = 0;
	
	for (unsigned int i=0; i<timeSeries_.size(); i++) {
		double max_i = getMaximumConcentration(timeSeries_.at(i));
		if (max_i > max)
			max = max_i;
	}
	return max;
}
	
// ----------------------------------------------------------------------------

/** Get the maximum concentration in the given time series. */
double TimeSeriesExperiment::getMaximumConcentration(const Mat_I_DP& ts) {

	double max = 0;
	
	for (int i=0; i<numTimePoints_; i++)
		for (int j=0; j<numGenes_; j++)
			if (ts[i][j] > max)
				max = ts[i][j];
	
	return max;
}

// ----------------------------------------------------------------------------

/**
 * Add experimental noise to the data
 */
void TimeSeriesExperiment::addNoise() {

	for (unsigned int i=0; i<timeSeries_.size(); i++)
		addNoise(timeSeries_.at(i));
	
	if (GnwSettings::Instance()->getModelTranslation())
		for (unsigned int i=0; i<timeSeriesProteins_.size(); i++)
			addNoise(timeSeriesProteins_.at(i));
	
	noiseHasBeenAdded_ = true;
}

// ----------------------------------------------------------------------------

/** Add experimental noise to the given data */
void TimeSeriesExperiment::addNoise(Mat_DP& ts) {
	
	for (int i=0; i<numTimePoints_; i++)
		for (int j=0; j<numGenes_; j++)
			ts[i][j] = Experiment::addNoise(ts[i][j]);
}

// ----------------------------------------------------------------------------

/** Normalize (i.e. divide by) the given maximum value */
void TimeSeriesExperiment::normalize(double max) {
	
	for (unsigned int i=0; i<timeSeries_.size(); i++)
		normalize(timeSeries_.at(i), max);
		
	if (GnwSettings::Instance()->getModelTranslation())
		for (unsigned int i=0; i<timeSeriesProteins_.size(); i++)
			normalize(timeSeriesProteins_.at(i), max);
}

// ----------------------------------------------------------------------------

/** Normalize (i.e. divide by) the given maximum value */
void TimeSeriesExperiment::normalize(Mat_DP& ts, double max) {
	
	for (int i=0; i<numTimePoints_; i++)
		for (int j=0; j<numGenes_; j++)
			ts[i][j] = ts[i][j]/max;
}
