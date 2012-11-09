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
*/

#ifndef RANDOMPARAMETERGAUSSIAN_H
#define RANDOMPARAMETERGAUSSIAN_H

#include "RandomParameter.h"

/** 
 * Generate random values from a Gaussian distribution within the interval
 * [min max].
 * 
 * If not set explicitly, the mean is the center of the interval (min+max)/2 and 
 * the standard deviation s (stdev_) so that the interval is 6 standard deviations
 * (max-min) = 6*s. This implies the following distribution of the generated values:
 * 
 *    ... mean     +1s     +2s     +3s (=max)
 *    ------|-------|-------|-------|--------->
 *    ...   | 34.1% | 13.6% | 2.1%  |
 *        
 * Set the appropriate flag to use a log-normal distributions (note, the  parameters
 * must be given on the linear scale, and the result is also returned on linear scale)
 * 
 * @author Daniel Marbach (firstname.name@gmail.com)
 * 
 */
class RandomParameterGaussian : public RandomParameter {

public:
	RandomParameterGaussian(double min, double max, bool logNormal);
	RandomParameterGaussian(double min, double max, double mean, double stdev, bool logNormal);
	~RandomParameterGaussian();
	
	double getRandomValue();

private:
	/** The mean */
	double mean_;
	/** The standard deviation */
	double stdev_;
	/** Use a log-normal distribution */
	bool logNormal_;
	
	void setMinMax(double min, double max);
	void setMeanStdev(double mean, double stdev);
	void setMean(double mean);
	
};

#endif