// $Id: GeneNetwork.h 29 2012-01-04 17:06:55Z jbao $

#ifndef GENENETWORK_H
#define GENENETWORK_H

#include <sstream>
#include <string>
#include <vector>
#include <sbml/SBMLTypes.h>
#include "nr/nr.h"
#include "Node.h"
#include "Edge.h"
#include "HillGene.h"
#include "GnwSettings.h"

/** This class represents a gene network.
 * 
 * It extends ImodNetwork, thus it contains all the
 * genes and the edges. Furthermore, it implements the state (gene and protein expression
 * levels) and functions to compute the production rates and to load and save the network
 * in SBML format.
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 * @author Daniel Marbach (firstname.name@gmail.com)
 * 
 */ 
class GeneNetwork { 

public:
	enum format {TSV, SBML};
	GeneNetwork();
	GeneNetwork(const GeneNetwork& grn);
 	GeneNetwork& operator=(const GeneNetwork& rhs);
	~GeneNetwork();

    static int gene_network_dynamics(double, const double xy[], double dxydt[], void *param);

	friend int gene_network_dynamics(double, const double xy[], double dxydt[], void *param);
	
	void load(const char *filename, format f);
	void setRank(const char *filename);
	std::vector<HillGene>& getRank();
	int getSize() { return nodes_.size(); }
	void writeTSV(const char *filename);
	void writeSBML(const char *filename);

	std::string getId() { return id_; }
	void setId(const std::string &id) { id_ = id; }
	
	int getIndexOfNode(std::string label);
	std::vector<HillGene>& getNodes() { return nodes_; }
	std::vector<Edge>& getEdges() { return edges_; }
	void getEdge(HillGene& src, HillGene& tgt, Edge& e);
	
	std::string getHeader(bool includeProteins);
	
	void computeDxydt(Vec_DP& xy, Vec_DP& dxydt);
	static void dynamics(Vec_I_DP& xy, Vec_O_DP& dxydt, GeneNetwork& grn);
	void computeInputs(Vec_I_DP& x, Vec_O_DP& productionRates);
	
	Vec_DP& getX() { return x_; }
	Vec_DP& getY() { return y_; }
	
	void randomInitialization();
	void prune(std::vector<HillGene>& genes);

    std::map<std::string, int> inDegrees;
    std::map<std::string, int> outDegrees;
    std::map<std::string, int> degrees;
    std::vector< std::pair<std::string, int> > inDegreesVec;
    std::vector< std::pair<std::string, int> > outDegreesVec;
    std::vector< std::pair<std::string, int> > degreesVec;

    void sanityCheck();

private:
	/** Current gene expression levels */
	Vec_DP x_;
	/** Current protein expression levels */
	Vec_DP y_;
	std::vector<HillGene> nodes_;
	std::vector<Edge> edges_;
	std::string id_;
	std::vector<HillGene> rank_;
	
	void load_tsv(const char *filename);
	void load_sbml(const char *filename);
	
	std::string getGeneReactantId(std::string& id);
	
	void computeMRnaProductionRates(Vec_DP& productionRates);
	
	void getInputs(HillGene& gene, std::vector<std::string>& inputs);
	void getOutputs(HillGene& gene, std::vector<std::string>& outputs);
	void getInputSigns(HillGene& gene, std::vector<Edge::sign>& inputSigns);
	void initializeInputWiring();
	
	void setEdgeTypesAccordingToDynamicalModel();
	
	int getIndexOfNode(HillGene& gene);
	
};



#endif
