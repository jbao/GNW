// $Id: main.cpp 29 2012-01-04 17:06:55Z jbao $

//#include <cstdlib>
//#include "logging/logging.h"
//using namespace ::logging;
#include <sstream>
#include <iostream>
#include "getopt_pp.h"
#include "GnwSettings.h"
#include "GeneNetwork.h"
#include "BenchmarkGenerator.h"

using namespace GetOpt;

int main(int argc, char **argv) {
	
    GetOpt_pp ops(argc, argv);
    std::string outdir, indir, type;
    int nid, size, perturb;
    double lower, upper;
    ops >> Option('o', "out-dir", outdir);
    ops >> Option('i', "in-dir", indir);
    ops >> Option('t', "type", type);
    ops >> Option('n', "network-id", nid);
    ops >> Option('s', "size", size);
    ops >> Option('p', "perturbation", perturb);

    //int i = 0;
    //char *task = getenv("SGE_TASK_ID");
    //if (task != NULL)
    //    i = atoi(task) - 1;
    //lower = i / 10.0;
    //upper = (i + 1) / 10.0;
    int i = perturb - 1;

    //std::string type = "ecoli";
	//std::string wd = "/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/gnw/prune/";
	//std::string wd = "/export/work/jbao/data/DREAM/gnw/" + type + 
    //    "/gnw/full/";
	//std::string wd = "/home/jbao/data/DREAM/gnw/scalefree/gnw/Size1000/";
	GnwSettings::Instance()->setPerturbationNumber(i + 1);
	GnwSettings::Instance()->setRandomSeed(42);
	GnwSettings::Instance()->setOutputDirectory(outdir);
	GnwSettings::Instance()->setModelTranslation(false);
	GnwSettings::Instance()->generateTsMultifactorial(true);
	GnwSettings::Instance()->generateTsConstantInput(false);
	GnwSettings::Instance()->setAddMicroarrayNoise(false);
	GnwSettings::Instance()->setNumTimeSeries(1);
	GnwSettings::Instance()->setMaxtTimeSeries(200);
	GnwSettings::Instance()->setDt(10);
	GnwSettings::Instance()->setMaxtSteadyStateODE(500);
	//GnwSettings::Instance()->setRelativePrecision(0);
	GnwSettings::Instance()->setPerturbationFraction(0.8, 1);
	
    GeneNetwork *net = new GeneNetwork();
	//net->load("/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/Networks/InSilicoSize100-Ecoli1.tsv", GeneNetwork::TSV);
    //net->load("/home/jbao/data/DREAM/DREAM3/prune/Size1000/InSilicoSize1000-Ecoli1.tsv", GeneNetwork::TSV);
	//net->load("/home/jbao/data/DREAM/gnw/ecoli/gnw/Size1000/default/ecoli-1_goldstandard_signed.tsv", GeneNetwork::TSV);
    std::stringstream ss, ssn, sssize;
    ss << i + 1;
    ssn << nid;
    sssize << size;
    std::string filename = indir + type + "-" + ssn.str() + "_" + sssize.str() + 
        ".tsv";
    //std::string filename = indir + "ecoli_transcriptional_network_regulonDB_6_2.tsv";
    //std::string filename = indir + "yeast_transcriptional_network_Balaji2006.tsv";
    //std::string filename = indir + "ecoli-full.xml";
    //std::string filename = indir + type + "-" + ss.str() + "_" + sssize.str() + 
    //    "_" + ssn.str() + ".xml";
    net->load(filename.c_str(), GeneNetwork::TSV);
    net->randomInitialization();
    net->sanityCheck();
	//net->load("/Users/bao/work/DREAM/gnw/ecoli/gnw/Size1000/default/rewiring/ecoli-1_1000_0.xml", GeneNetwork::SBML);
	//net->load("/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/gnw/prune/full.xml", GeneNetwork::SBML);
	//net->setRank("/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/gnw/mds/ranked_ecoli1_100");
	//net->setRank("/home/jbao/data/DREAM/DREAM3/prune/ranked_ecoli1_100");
	//std::vector<HillGene> all_nodes = net->getNodes();
	//for (int i = 0; i < res.size(); ++i) {
	//	GeneNetwork n = GeneNetwork(*net);
		//std::vector<HillGene> toPrune;
		//toPrune.push_back(all_nodes.at(i));
		//net->prune(toPrune);
		//std::stringstream ss;
		//ss << i + 1;
	//net->setId(type + "-" + ss.str() + "_" + ssn.str());
	net->setId(type + "-" + ssn.str());
    //net->setId(type + "-full");

    //for (std::vector< std::pair<std::string,int> >::iterator it = net->outDegreesVec.begin();
    //     it != net->outDegreesVec.end(); ++it)
	//    std::cerr << it->first << " " << it->second << std::endl;
	
	GnwSettings::Instance()->setRandomSeed(i);
	//for (int ii = 0; ii < i; ++ii)
    //    int tmp = GnwSettings::Instance()->getNormalDistributionNumber(0, 1);
    BenchmarkGenerator *bg = new BenchmarkGenerator(*net);
    bg->generateGoldStandard();
    delete bg;
	//}
	
	delete net;
	//delete bg;
	
	return 0;
}
