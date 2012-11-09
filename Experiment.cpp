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

$Id: Experiment.cpp 19 2011-06-17 08:20:59Z jbao $
*/

#include <cmath>
#include "GnwSettings.h"
#include "Exception.h"
#include "Experiment.h"	
	
// ============================================================================
// PUBLIC METHODS

/**
 * Default constructor
 */
Experiment::Experiment(const std::string& label) { 
	GnwSettings *set = GnwSettings::Instance();
	//solverType_ = solverType;
	label_ = label;
	//grn_ = null;
	numGenes_ = 0;
	perturbation_ = new PerturbationMultifactorial();
	
	numExperiments_ = 1;
	//xy0_ = null;
	//modelTranslation_ = set.getModelTranslation();
	//normalDistribution_ = set.getNormalDistribution();
	
	addNormalNoise_ = set->getAddNormalNoise();
	addLognormalNoise_ = set->getAddLognormalNoise();
	addMicroarrayNoise_ = set->getAddMicroarrayNoise();
	if ((addMicroarrayNoise_ && addNormalNoise_) || (addMicroarrayNoise_ && addLognormalNoise_))
		throw Exception("You can't add both normal/lognormal noise and microarray noise");

	normalStdev_ = set->getNormalStdev();
	lognormalStdev_ = set->getLognormalStdev();
	
	noiseHasBeenAdded_ = false;
}
	
// ---------------------------------------------------------------------------
 
Experiment::Experiment(Perturbation *perturbation, std::string label) {
	
	GnwSettings *set = GnwSettings::Instance();
	//solverType_ = solverType;
	label_ = label;
	//grn_ = null;
	numGenes_ = 0;
	perturbation_ = perturbation;
	if (perturbation_->getNumGenes() != 0)
		numExperiments_ = perturbation_->getNumPerturbations();
	else
		numExperiments_ = 1;
	//xy0_ = null;
	//modelTranslation_ = set.getModelTranslation();
	//normalDistribution_ = set.getNormalDistribution();
	
	addNormalNoise_ = set->getAddNormalNoise();
	addLognormalNoise_ = set->getAddLognormalNoise();
	addMicroarrayNoise_ = set->getAddMicroarrayNoise();
	if ((addMicroarrayNoise_ && addNormalNoise_) || (addMicroarrayNoise_ && addLognormalNoise_))
		throw Exception("You can't add both normal/lognormal noise and microarray noise");

	normalStdev_ = set->getNormalStdev();
	lognormalStdev_ = set->getLognormalStdev();
	
	noiseHasBeenAdded_ = false;
	
}
	
// ----------------------------------------------------------------------------

Experiment::Experiment(const Experiment& e) {
	label_ = e.label_;
	numGenes_ = e.numGenes_;
	perturbation_ = e.perturbation_;
	numExperiments_ = e.numExperiments_;
	addNormalNoise_ = e.addNormalNoise_;
	addLognormalNoise_ = e.addLognormalNoise_;
	addMicroarrayNoise_ = e.addMicroarrayNoise_;
	normalStdev_ = e.normalStdev_;
	lognormalStdev_ = e.lognormalStdev_;
	noiseHasBeenAdded_ = e.noiseHasBeenAdded_;
    grn_ = e.grn_;
}

// ----------------------------------------------------------------------------

Experiment& Experiment::operator=(const Experiment& rhs) {
	label_ = rhs.label_;
	numGenes_ = rhs.numGenes_;
	perturbation_ = rhs.perturbation_;
	numExperiments_ = rhs.numExperiments_;
	addNormalNoise_ = rhs.addNormalNoise_;
	addLognormalNoise_ = rhs.addLognormalNoise_;
	addMicroarrayNoise_ = rhs.addMicroarrayNoise_;
	normalStdev_ = rhs.normalStdev_;
	lognormalStdev_ = rhs.lognormalStdev_;
	noiseHasBeenAdded_ = rhs.noiseHasBeenAdded_;
	grn_ = rhs.grn_;

	return *this;
}
	
// ----------------------------------------------------------------------------
	
Experiment::~Experiment() {

}	

// ----------------------------------------------------------------------------
	
/** Add log-normal noise to the data point x, set values below threshold to zero. TODO check that normal/lognormal works */
double Experiment::addNoise(double x) {
	
	if (x < 0)
		throw Exception("Experiment:addNoise(): x < 0!");
	
	double xPlusNoise = x;
	
	// Note, in the constructor we tested that not both microarray noise and normal/lognormal noise
	// is added, which makes no sense. However, normal and lognormal noise can be added together
	if (addLognormalNoise_)
		xPlusNoise = addLogNormalNoise(x);
	if (addNormalNoise_)
		xPlusNoise += GnwSettings::Instance()->getRNG()->randNormInv(0, normalStdev_);
	if (addMicroarrayNoise_)		
		xPlusNoise = addMicroarrayNoise(x);
		
	if (xPlusNoise < 0)
		xPlusNoise = 0;
	
	return xPlusNoise;
}

// ----------------------------------------------------------------------------
	
/** Add log-normal noise to the data point x */
double Experiment::addLogNormalNoise(double x) {
	
	if (x < 0)
		throw Exception("Experiment:addLogNormalNoise(): x < 0!");
	else if (x == 0.0)
		return 0.0;
	
	// transform to log-scale
	double xLog = log10(x); // TODO does it make a difference whether we use log10, log2...?
	
	// add normal noise with mean zero: y = x + n(0, s) = n(m, s)
	xLog = GnwSettings::Instance()->getRNG()->randNormInv(xLog, lognormalStdev_);
	
	return pow(10., xLog);
}


//	----------------------------------------------------------------------------

/** Use the model of microarray noise by Tu et al. (2002) */
double Experiment::addMicroarrayNoise(double x) {
	
	if (x < 0)
		throw Exception("Experiment:addMicroarrayNoise(): x < 0!");
	else if (x == 0.0)
		return 0.0;
	
	// TODO allow these parameters to be set by the user? (btw, this could be done only once)
	double alpha = 0.001;
	double beta = 0.69;
	double K = 0.01;
	
	double variance = alpha + (beta - alpha)/(1 + (x/K));
	double w = GnwSettings::Instance()->getRNG()->randNormInv(0, sqrt(variance));
	
	return x*exp(w);
}

// ---------------------------------------------------------------------------

void Experiment::rkdumb(Vec_I_DP &vstart, const DP x1, const DP x2,
	void (GeneNetwork::*derivs)(const DP, Vec_I_DP &, Vec_O_DP &)) {

	int i,k;
	DP x,h;

	Vec_DP &xx=*xx_p;
	Mat_DP &y=*y_p;
	int nvar=y.nrows();
	int nstep=y.ncols()-1;
	Vec_DP v(nvar),vout(nvar),dv(nvar);
	for (i=0;i<nvar;i++) {
		v[i]=vstart[i];
		y[i][0]=v[i];
	}
	xx[0]=x1;
	x=x1;
	h=(x2-x1)/nstep;
	for (k=0;k<nstep;k++) {
		(grn_.*derivs)(x,v,dv);
		rk4(v,dv,x,h,vout,derivs);
		if (x+h == x)
			NR::nrerror("Step size too small in routine rkdumb");
		x += h;
		xx[k+1]=x;
		for (i=0;i<nvar;i++) {
			v[i]=vout[i];
			y[i][k+1]=v[i];
		}
	}
}

// ---------------------------------------------------------------------------

void Experiment::rk4(Vec_I_DP &y, Vec_I_DP &dydx, const DP x, const DP h,
	Vec_O_DP &yout, void (GeneNetwork::*derivs)(const DP, Vec_I_DP &, Vec_O_DP &)) {
	
	int i;
	DP xh,hh,h6;

	int n=y.size();
	Vec_DP dym(n),dyt(n),yt(n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(grn_.*derivs)(xh,yt,dyt);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(grn_.*derivs)(xh,yt,dym);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(grn_.*derivs)(x+h,yt,dyt);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

