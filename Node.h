#ifndef NODE_H
#define NODE_H

#include <string>
#include <map>
#include <vector>
#include "RegulatoryModule.h"
//#include "GeneNetwork.h"

class Node {

public:
	Node() {};
	Node(std::string label) : label_(label) {};
	virtual ~Node() {};
	inline bool operator== (const Node &rhs) const {return label_ == rhs.label_;}
	
	std::string getLabel() {return label_;}
	void setLabel(std::string label) {label_ = label;}
	
	//virtual void initialization(std::map<std::string, double> params, std::vector<HillGene> inputGenes) {};
	
	virtual void compileParameters(std::map<std::string, double> params) {}
	
	//virtual std::vector<Node> getInputGenes() { return inputGenes_; }
 	//virtual void setInputGenes(std::vector<Node> inputs) {}
 	
 	//virtual std::vector<RegulatoryModule> getRegulatoryModules() { return regulatoryModules_; }
 	//virtual void setRegulatoryModules(std::vector<RegulatoryModule> modules) {}
	
protected:
	//std::vector<Node> inputGenes_;
	//std::vector<RegulatoryModule> regulatoryModules_;
	std::string label_;
	//GeneNetwork grn_;

};

#endif