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

$Id: SteadyStateExperiment.cpp 27 2011-09-19 08:37:15Z jbao $
*/

#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include "GnwSettings.h"
#include "SteadyStateExperiment.h"
#include "logging/logging.h"
using namespace ::logging;
	
// ---------------------------------------------------------------------------
	
SteadyStateExperiment::SteadyStateExperiment(const std::string& label) : Experiment(label) {
 	maxtODE_ = GnwSettings::Instance()->getMaxtSteadyStateODE();
 	//ssPerturbation_ = Mat_DP(numExperiments_, numGenes_);
 	
 	std::stringstream ss;
 	ss << maxtODE_;
 	//::logging::log::emit<Debug>() << "maxt = " << ss.str().c_str() << ::logging::log::endl;
}
 	
// ---------------------------------------------------------------------------
 	
SteadyStateExperiment::SteadyStateExperiment(PerturbationMultifactorial *perturbation, std::string label) : Experiment(perturbation, label) {
 	maxtODE_ = GnwSettings::Instance()->getMaxtSteadyStateODE();
 	//ssPerturbation_ = Mat_DP(numExperiments_, numGenes_);
}
	
// ----------------------------------------------------------------------------

/**
 * Run all steady-state experiments
 */
void SteadyStateExperiment::run(Vec_DP& xy0) {
	
	xy0_ = xy0;
	
	try {
		std::string simulationType = "ODEs";
		//if (solverType_ == Solver.type.SDE)
		//	simulationType = "SDEs";
		::logging::log::emit<Info>() << "Simulating steady-state " << label_.c_str() << " using " << simulationType.c_str() << " ..." << ::logging::log::endl;
		
		//if (solverType_ == Solver.type.SDE && timeToConvergenceODE_ == null && maxtSDE_ < 0)
		//	throw new RuntimeException("For SDE steady-state simulation, either specify timeToConvergenceODE_ or maxtSDE_");
			
		Mat_DP ssPerturbation_(numExperiments_, numGenes_);
		if (GnwSettings::Instance()->getModelTranslation())
			Mat_DP ssPerturbationProteins_(numExperiments_, numGenes_);
			
		//if (solverType_ == Solver.type.ODE)
		//	timeToConvergenceODE_ = new ArrayList<Double>();
		
		computeSteadyStates();
		
		// display the longest time to convergence
		if (simulationType == "ODEs") {
			double max = -1;
			for (unsigned int i=0; i<timeToConvergenceODE_.size(); i++)
				if (timeToConvergenceODE_[i] > max)
					max = timeToConvergenceODE_[i];
			std::stringstream ss;
			ss << max;
			::logging::log::emit<Info>() << "Duration of the longest steady state experiment = " << ss.str().c_str() << ::logging::log::endl;
		}
		//log.log(Level.INFO, ""); // empty line

	} catch (std::exception& e) {
		::logging::log::emit<Info>() << "SteadyStateExperiment::runAll(): Exception " << e.what() << ::logging::log::endl;
		//throw new RuntimeException();
	} 
}

// ----------------------------------------------------------------------------

/**
 * Run all steady-state experiments with 
 */
void SteadyStateExperiment::run() {
	
	//xy0_ = xy0;
	
	try {
		std::string simulationType = "ODEs";
		//if (solverType_ == Solver.type.SDE)
		//	simulationType = "SDEs";
		::logging::log::emit<Info>() << "Simulating steady-state " << label_.c_str() << " using " << simulationType.c_str() << " ..." << ::logging::log::endl;
		
		//if (solverType_ == Solver.type.SDE && timeToConvergenceODE_ == null && maxtSDE_ < 0)
		//	throw new RuntimeException("For SDE steady-state simulation, either specify timeToConvergenceODE_ or maxtSDE_");
			
		Mat_DP ssPerturbation_(numExperiments_, numGenes_);
		if (GnwSettings::Instance()->getModelTranslation())
			Mat_DP ssPerturbationProteins_(numExperiments_, numGenes_);
			
		//if (solverType_ == Solver.type.ODE)
		//	timeToConvergenceODE_ = new ArrayList<Double>();
		
		computeSteadyStates();
		
		// display the longest time to convergence
		if (simulationType == "ODEs") {
			double max = -1;
			for (unsigned int i=0; i<timeToConvergenceODE_.size(); i++)
				if (timeToConvergenceODE_[i] > max)
					max = timeToConvergenceODE_[i];
			std::stringstream ss;
			ss << max;
			::logging::log::emit<Info>() << "Duration of the longest steady state experiment = " << ss.str().c_str() << ::logging::log::endl;
		}
		//log.log(Level.INFO, ""); // empty line

	} catch (Exception& e) {
		::logging::log::emit<Error>() << "SteadyStateExperiment::run(): " << e.what().c_str() << ::logging::log::endl;
		//throw new RuntimeException();
	} 
}

// ----------------------------------------------------------------------------

/**
 * Print the mRNA data, append the given string to the filenames (e.g. "_nonoise_wildtype").
 * If translation is modelled, protein data is also printed
 */
void SteadyStateExperiment::printAll(std::string postfix) {
	printMRNA(postfix);
	if (GnwSettings::Instance()->getModelTranslation())
		printProteins(postfix);
}


// ----------------------------------------------------------------------------

/**
 * Print the mRNA data, append the given string to the filenames (e.g. "-nonoise"). 
 */
void SteadyStateExperiment::printMRNA(std::string postfix) {
	
	std::string absPath = GnwSettings::Instance()->getOutputDirectory();
	
	if (label_ != "wildtype" && postfix != "_noexpnoise")
		printSteadyStates(absPath + grn_.getId() + postfix + "_" + label_ + ".tsv", ssPerturbation_);
}


// ----------------------------------------------------------------------------

/**
 * Print the protein data, append the given string to the filenames (e.g. "-nonoise"). 
 */
void SteadyStateExperiment::printProteins(std::string postfix) {
	
	std::string absPath = GnwSettings::Instance()->getOutputDirectory();
	
	//if (!GnwSettings::Instance()->getModelTranslation())
	//	throw new IllegalArgumentException("SteadyStateExperiment:printProteins(): protein translation was not modeled");

	printSteadyStates(absPath + grn_.getId() + postfix + "_proteins_" + label_ + ".tsv", ssPerturbationProteins_);
}


// ----------------------------------------------------------------------------

/** 
 * Compute the steady-states for all the single-gene perturbations, one after the other.
 * The result is stored in ssPerturbation.  
 * @throws Exception 
 */ 
void SteadyStateExperiment::computeSteadyStates() {
	
	// apply each perturbation, one after the other, and compute the steady-states
	for (int i=0; i<numExperiments_; i++) {
		
		// the time limit for the simulation
		double maxt;
		//if (solverType_ == Solver.type.SDE) {
		//	if (maxtSDE_ > 0)
		//		maxt = maxtSDE_;
		//	else
		//		maxt = timeToConvergenceODE_.get(i); 
		//} else
		maxt = maxtODE_;
		
		// apply the perturbation
		//if (perturbation_->getNumGenes() != 0)
		//	perturbation_->applyPerturbation(i);
		
		// compute the steady-state
		computeSteadyState(maxt);
		
		// remove the perturbation
		//if (perturbation_->getNumGenes() != 0)
		//	perturbation_->restoreWildType();
		
		// put the steady-state into the corresponding line in ssPerturbation_
		ssPerturbation_ = Mat_DP(numExperiments_, numGenes_);
		Vec_DP x = grn_.getX();
		for (int j=0; j<numGenes_; j++)
			ssPerturbation_[i][j] = x[j];
		
		if (GnwSettings::Instance()->getModelTranslation()) {
			Vec_DP y = grn_.getY();
			for (int j=0; j<numGenes_; j++)
				ssPerturbationProteins_[i][j] = y[j];
		}
	}
	// remove the perturbation from the network
	//if (perturbation_->getNumGenes() != 0)
	//	perturbation_->restoreWildType();
}


// ============================================================================
// PRIVATE METHODS

/**
 * Compute the steady state of the network after integrating from the given
 * initial conditions x0 and y0.
 * @throws Exception 
 */
void SteadyStateExperiment::computeSteadyState(double maxt) {
					
	Vec_DP xy0;
	constructInitialCondition(xy0); // initial condition
	t_ = 0;
	double dt = GnwSettings::Instance()->getDt();
	
	//if (dt <= 0 || dt > maxt)
	//	throw new IllegalArgumentException("Interval between two measuread points must be >0 and <maxt.");
	
	// TODO: verify why GnwSettings->dt_ is used => most likely to be removed
	// replace that by being sure that dt is a multiple of the internal integration step-size
	// as dt is calculated from maxt_ and numTimePoints there must be no problem (see reported bug)
	double frac = maxt/dt;
	if (frac - static_cast<int>(frac) != 0) 
		throw Exception("Duration (maxt) must be a multiple of numTimePoints-1.");
	
	//Solver solver = new Solver(solverType_, grn_, xy0);
	Vec_DP previousState = xy0;
	Vec_DP state = xy0;
	do {
		std::stringstream ss;
 		ss << t_;
 		//::logging::log::emit<Debug>() << "t = " << ss.str().c_str() << 
        //    ::logging::log::endl;
	
		//double t1 = t_;
		// this steps the time by dt_, but using a smaller internal step size of the solver
		// (getRate() may be called several times for one step)
		t_ += dt;
		Vec_DP dxydt(numGenes_);
		grn_.computeDxydt(previousState, dxydt);
		for (int i = 0; i < numGenes_; ++i) {
			state[i] = previousState[i] + dxydt[i]*dt;
            if (state[i] < 0)
                state[i] = 0;
        }
		
		//if (t_ != t1 + dt)
		//	throw new RuntimeException("Solver failed to step time by dt, expected t = " + (t1+dt) + ", obtained t = " + t_);
		
	} while (!converged(previousState, state) && t_ < maxt);
	// note, the state at the last step is already saved both in ODE.state and grn.x_, grn.y_
	
	// save the time of this experiment
	//if (solverType_ == Solver.type.ODE)
	timeToConvergenceODE_.push_back(t_);
	
	// Check the max rate of change at the found solution
	Vec_DP lastX = grn_.getX();
	Vec_DP lastY = grn_.getY();
	Vec_DP dxydt(numGenes_);
	Vec_DP xy(numGenes_);
	if (GnwSettings::Instance()->getModelTranslation()) {
		Vec_DP dxydt(2*numGenes_);
		Vec_DP xy(lastX.size()+lastY.size(), 0);
		for (int i = 0; i < lastX.size(); ++i)
			xy[i] = lastX[i];
		for (int i = lastX.size(); i < xy.size(); ++i)
			xy[i] = lastY[i];
	} else {
		//Vec_DP xy(numGenes_);
		//Vec_DP dxydt(numGenes_);
		xy = lastX;
		//std::stringstream ss;
 		//ss << xy[0];
 		//::logging::log::emit<Debug>() << "xy[0] = " << ss.str().c_str() << ::logging::log::endl;
	}
		
	grn_.computeDxydt(xy, dxydt);
	
	double max = 0;
	for (int i=0; i<dxydt.size(); i++)
		if (dxydt[i] > max)
			max = dxydt[i];

	std::stringstream ss_t, ss_max;
	ss_t << t_;
	ss_max << max;
	::logging::log::emit<Info>() << "Saved state at t = " << ss_t.str().c_str() << ", with maximum dx_i/dt = " << ss_max.str().c_str() << ::logging::log::endl;
	
    steadyStates_ = xy;
	//getSDESolver()
	//if (solverType_ == Solver.type.SDE && solver.getXNegativeCounter() > 0)
	//	log.log(Level.INFO, "SDE: " + solver.getXNegativeCounter() + " times a concentration became negative due to noise and was set to 0");
}


// ----------------------------------------------------------------------------

/**
 * Construct xy0, the initial conditions as an array of double. If x0 is null,
 * we estimate the initial conditions as the concentration of the genes in the
 * absence of regulation.
 */
void SteadyStateExperiment::constructInitialCondition(Vec_DP& xy0) {
	
	//Vec_DP xy0;
	
	if (xy0_.size() == 0) {
		// construct a vector of zeros
		Vec_DP zeros;
		if (GnwSettings::Instance()->getModelTranslation()) {
			xy0 = Vec_DP(2*numGenes_);
			zeros = Vec_DP(0.0, 2*numGenes_);
			for (int i = 0; i < 2*numGenes_; ++i)	
                zeros[i] = GnwSettings::Instance()->getUniformDistributionNumber(1);
		} else {
			xy0 = Vec_DP(numGenes_);
			zeros = Vec_DP(0.0, numGenes_);
			for (int i = 0; i < numGenes_; ++i)	
                zeros[i] = GnwSettings::Instance()->getUniformDistributionNumber(1);
		}
		//for (int i=0; i<zeros.size(); i++)
		//	zeros[i] = 0;

		// Estimate the initial conditions as the concentration of the genes without regulation.
		// 0 = m*f(0) -delta*x_i  =>  x_i = m*f(0) / delta
		grn_.computeDxydt(zeros, xy0); // since x=0, this actually computes m*f(0)
		
		for (int i=0; i<numGenes_; i++) {
			xy0[i] /= grn_.getNodes().at(i).getDelta();
            if (xy0[i] < 0)
                xy0[i] = 0;
		}

		// Corresponding initial conditions for the proteins:
		// 0 = mTranslation*x_i - deltaProt*y_i  =>  y_i = mTranslation*x_i / deltaProt
		if (GnwSettings::Instance()->getModelTranslation()) {
			for (int i=0; i<numGenes_; i++) {
				double m = grn_.getNodes().at(i).getMaxTranslation();
				double d = grn_.getNodes().at(i).getDeltaProtein();
				xy0[numGenes_+i] = m*xy0[i] / d;
			}
		}
		//xy0_ = new DenseDoubleMatrix1D(xy0.length);
		//xy0_.assign(xy0);
	
	} else {
		xy0 = xy0_;
	}
	//return xy0;
}

// ----------------------------------------------------------------------------

void SteadyStateExperiment::fdjac(Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df,
    void (*dynamics)(Vec_I_DP&, Vec_O_DP&, GeneNetwork&)) {
    
    const DP EPS=1.0e-8;
    int i,j;
    DP h,temp;

    int n=x.size();
    Vec_DP f(n);
    for (j=0;j<n;j++) {
        temp=x[j];
        h=EPS*fabs(temp);
        if (h == 0.0) h=EPS;
        x[j]=temp+h;
        h=x[j]-temp;
        dynamics(x,f, grn_);
        x[j]=temp;
        for (i=0;i<n;i++)
            df[i][j]=(f[i]-fvec[i])/h;
    }
}

// ----------------------------------------------------------------------------

/**
 * Using analytical solutions in Jacobian.
 */
void SteadyStateExperiment::jac(Vec_I_DP &xvec, Mat_O_DP &df) {
	::logging::log::emit<Info>() << "Calculating analytical Jacobian matrix at the steady state..." << ::logging::log::endl;
    
    // diagonal elements of the Jacobian are the degradation rates
    assert(numGenes_ == df.nrows() && numGenes_ == df.ncols());
    for (int i = 0; i < numGenes_; ++i) {
        df[i][i] = - grn_.getNodes().at(i).getDelta();
        //std::stringstream iss,dfss;
        //iss << i;
        //dfss << df[i][i];
	    //::logging::log::emit<Debug>() << "jac[" << iss.str().c_str() << "][" << 
        //    iss.str().c_str() << "] = " << dfss.str().c_str() << ::logging::log::endl;
    }

    // entries for input genes are defined analytically, otherwise 0
    for (int i = 0; i < numGenes_; ++i) {
        std::vector<std::string> inputGenes = grn_.getNodes().at(i).getInputGenes();
        int numInputs = grn_.getNodes().at(i).getNumInputs();
        std::map<std::string,double> params;
        grn_.getNodes().at(i).compileParameters(params);
	    int numStates = static_cast<int>(pow(2.0, numInputs));
        
        for (int j = 0; j < numInputs; ++j) {
            //std::stringstream jss;
            //jss << j;
            //::logging::log::emit<Debug>() << "\tInput " << jss.str().c_str() << ::logging::log::endl;
            double sum = 0;
            int iInput = grn_.getIndexOfNode(inputGenes[j]);
		    
            for (int ns = 0; ns < numStates; ++ns) {
                //std::stringstream nss;
                //nss << ns;
                //::logging::log::emit<Debug>() << "\t\tState " << nss.str().c_str() << ::logging::log::endl;
                char bits[numInputs];
                grn_.getNodes().at(i).dec2bin(ns, bits);
                std::stringstream ss;
                ss << bits;
                std::string s = ss.str();
		        double p = 1; // the probability of being in state ns
                std::stringstream ass;
                ass << "a_" << ns;
                double alpha = params[ass.str()];

                for (int ni = 0; ni < numInputs; ++ni) {
                    //std::stringstream niss;
                    //niss << ni;
                    //::logging::log::emit<Debug>() << "\t\t\tInput " << niss.str().c_str() << ::logging::log::endl;
                    int idx = grn_.getIndexOfNode(inputGenes[ni]);
                    double x = xvec[idx];
                    std::stringstream kss,nss;
                    kss << "k_" << ni + 1;
                    double k = params[kss.str()];
                    nss << "n_" << ni + 1;
                    double n = params[nss.str()];

			        // p = alpha * (1 - H) for '1'/binding, p = alpha * H otherwise
                    // s.length() has to be cast to int
                    if (static_cast<int>(s.length())-ni-1 >= 0 && 
                            s.at(static_cast<int>(s.length())-ni-1) == '1') {
                        if (ni == j) { // calculate the derivative
                            if (x == 0) 
                                p *= 0;
                            else
                                p *= n * pow(x/k,n) / (x * pow(pow(x/k,n) + 1, 2));
                        }
                        else // keep the term
                            p *= 1 - 1 / (pow(x/k,n) + 1);
                    }
			        else {
                        if (ni == j) { // calculate the derivative
                            if (x == 0)
                                p *= 0;
                            else
                                p *= -n * pow(x/k,n) / (x * pow(pow(x/k,n) + 1, 2));
                        }
                        else // keep the term
                            p *= 1 / (pow(x/k,n) + 1);
                    }
                    sum += alpha * p;
                }
                //assert(p >= 0 && p <= 1);
                //assert((sum += p) >= 0); // always true, just to compute the sum
            }
            df[iInput][i] = sum * grn_.getNodes().at(i).getMaxTranscription();
            //std::stringstream iss,iiss,dfss;
            //iss << i;
            //iiss << iInput;
            //dfss << df[iInput][i];
            //::logging::log::emit<Debug>() << "jac[" << iiss.str().c_str() << "][" << 
            //    iss.str().c_str() << "] = " << dfss.str().c_str() << ::logging::log::endl;
        }
    }
}

// ----------------------------------------------------------------------------

/**
 * Print the given steady-state vector (data). This function is usually used to print the wild-type.
 */
void SteadyStateExperiment::printSteadyStates(std::string filename, const Vec_DP& wt) {
	
	// copy wt to a 2D matrix
	Mat_DP data(1, wt.size());
	for (int i=0; i<wt.size(); i++)
		data[0][i] = wt[i];
	
	printSteadyStates(filename, data);
}


// ----------------------------------------------------------------------------

/**
 * Print the perturbation experiments (data).
 */
void SteadyStateExperiment::printSteadyStates(std::string filename, const Mat_DP& data) {
	
	try {
		::logging::log::emit<Info>() << "Writing file " << filename.c_str() << ::logging::log::endl;
		ofstream data_file(filename.c_str()); 
		std::stringstream ss;

		data_file << grn_.getHeader(false);
		
		// Data
		if (&data != NULL) {
			for (int i=0; i<data.nrows(); i++) {
				//fw.write("\"" + grn_.getNode(i).getLabel() + type + "\"");

				for (int j=0; j<data.ncols()-1; j++) {
					ss.str("");
					ss << data[i][j];
					data_file << setprecision(7) << ss.str() << "\t";//Double.toString(value));
				}
				ss.str("");
				ss << data[i][data.ncols()-1];
				data_file << setprecision(7) << ss.str() << std::endl;
			}
		}

		// Close file
		data_file.close();

	} catch (std::exception& e) {
		::logging::log::emit<Info>() << "SteadyStateExperiment::printSteadyStates(): " << e.what() << ::logging::log::endl;
		//throw new RuntimeException();
	}
}

// ----------------------------------------------------------------------------

/**
 * Print the Jacobian matrix.
 */
void SteadyStateExperiment::printJacobian(std::string filename, const Mat_DP& data) {
	
	try {
		::logging::log::emit<Info>() << "Writing file " << filename.c_str() << ::logging::log::endl;
		ofstream data_file(filename.c_str()); 
		std::stringstream ss;

        // header
		data_file << grn_.getHeader(false);

		// matrix
		if (&data != NULL) {
			for (int i=0; i<data.nrows(); i++) {
				//fw.write("\"" + grn_.getNode(i).getLabel() + type + "\"");

				for (int j=0; j<data.ncols()-1; j++) {
					ss.str("");
					ss << data[i][j];
					data_file << setprecision(7) << ss.str() << "\t";//Double.toString(value));
				}
				ss.str("");
				ss << data[i][data.ncols()-1];
				data_file << setprecision(7) << ss.str() << std::endl;
			}
		}

		// Close file
		data_file.close();

	} catch (std::exception& e) {
		::logging::log::emit<Info>() << "SteadyStateExperiment::printJacobian(): " << e.what() << ::logging::log::endl;
		//throw new RuntimeException();
	}
}

// ----------------------------------------------------------------------------

/** 
 * Implements the method gsl_multiroot_test_delta() of GSL:
 * This function tests for the convergence of the sequence by comparing the last step dx
 * with the absolute error epsabs and relative error epsrel to the current position x.
 * The test returns true if the following condition is achieved:
 * 		|dx_i| < epsabs + epsrel |x_i|
 * for each component of x and returns false otherwise.
 */
bool SteadyStateExperiment::converged(Vec_DP& previousState, Vec_DP& state) {
	 
	double absolutePrecision = GnwSettings::Instance()->getAbsolutePrecision();
	double relativePrecision = GnwSettings::Instance()->getRelativePrecision();
	 
	for (int i=0; i<previousState.size(); i++) {
		
		double dxy = abs(previousState[i] - state[i]); 
		
		if (dxy > absolutePrecision + relativePrecision*abs(state[i])) {
			// remember point
			for (int k=0; k<previousState.size(); k++)
				previousState[k] = state[k];
			
			return false;
		}
	}
	return true;
}
	
// ----------------------------------------------------------------------------

/**
 * Get the maximum concentration of all experiments.
 * For now, only get the max concentration between mRNA levels. Later, perhaps
 * leave the choice to the user if he prefer, e.g., to normalise by the mRNA
 * OR¬†protein levels.
 */
double SteadyStateExperiment::getMaximumConcentration() {

	double max = 0;
	for (int i=0; i<numExperiments_; i++) {
		for (int j=0; j<numGenes_; j++) {
			if (ssPerturbation_[i][j] > max)
				max = ssPerturbation_[i][j];
		}
	}
	
	return max;
}

// ----------------------------------------------------------------------------

/**
 * Add log-normal noise to the data 
 */
void SteadyStateExperiment::addNoise() {
	
	GnwSettings *settings = GnwSettings::Instance();
	
	// mRNA
	for (int i=0; i<numExperiments_; i++)
		for (int j=0; j<numGenes_; j++)
			ssPerturbation_[i][j] = Experiment::addNoise(ssPerturbation_[i][j]);
	
	// proteins
	if (settings->getModelTranslation()) {
		for (int i=0; i<numExperiments_; i++)
			for (int j=0; j<numGenes_; j++)
				ssPerturbationProteins_[i][j] = Experiment::addNoise(ssPerturbationProteins_[i][j]);
	}
	
	noiseHasBeenAdded_ = true;
}

// ----------------------------------------------------------------------------

/**
 * Normalize (i.e. divide by) the given maximum value 
 */
void SteadyStateExperiment::normalize(double max) {
	
	GnwSettings *settings = GnwSettings::Instance();

	// mRNA
	for (int i=0; i<numExperiments_; i++)
		for (int j=0; j<numGenes_; j++)
			ssPerturbation_[i][j] = ssPerturbation_[i][j]/max;

	// proteins
	if (settings->getModelTranslation()) {
		for (int i=0; i<numExperiments_; i++)
			for (int j=0; j<numGenes_; j++)
				ssPerturbationProteins_[i][j] = ssPerturbationProteins_[i][j]/max;
	}
	
}
