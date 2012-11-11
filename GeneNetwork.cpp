// $Id: GeneNetwork.cpp 29 2012-01-04 17:06:55Z jbao $

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include "Exception.h"
#include "GeneNetwork.h"
#include "HillGene.h"
#include "RegulatoryModule.h"
#include "logging/logging.h"
using namespace ::logging;
		
LIBSBML_CPP_NAMESPACE_USE
	
struct IntCmp {
    bool operator()(const std::pair<std::string,int>& lhs, 
                    const std::pair<std::string,int>& rhs) {
        return lhs.second < rhs.second;
    }
};

// ============================================================================
// INITIALIZATION

/**
 * Default constructor
 */
GeneNetwork::GeneNetwork() {
		//x_ = NULL;
		//y_ = NULL;
		id_ = "grn";
}

// ----------------------------------------------------------------------------

GeneNetwork::GeneNetwork(const GeneNetwork& grn) {
	x_ = grn.x_;
	y_ = grn.y_;
	nodes_ = grn.nodes_;
	edges_ = grn.edges_;
	id_ = grn.id_;
    inDegrees = grn.inDegrees;
    outDegrees = grn.outDegrees;
    degrees = grn.degrees;
    inDegreesVec = grn.inDegreesVec;
    outDegreesVec = grn.outDegreesVec;
    degreesVec = grn.degreesVec;
}

// ----------------------------------------------------------------------------

GeneNetwork& GeneNetwork::operator=(const GeneNetwork& rhs) {
	x_ = rhs.x_;
	y_ = rhs.y_;
	nodes_ = rhs.nodes_;
	edges_ = rhs.edges_;
	id_ = rhs.id_;
    inDegrees = rhs.inDegrees;
    outDegrees = rhs.outDegrees;
    degrees = rhs.degrees;
    inDegreesVec = rhs.inDegreesVec;
    outDegreesVec = rhs.outDegreesVec;
    degreesVec = rhs.degreesVec;
	
	return *this;
}

// ---------------------------------------------------------------------------

GeneNetwork::~GeneNetwork() {
	//for (std::vector<HillGene>::iterator it = nodes_.begin(); it != nodes_.end(); ++it)
	//	delete *it;
		
	//for (std::vector<Edge>::iterator it = edges_.begin(); it != edges_.end(); ++it)
	//	delete *it;
		
	nodes_.clear();
	edges_.clear();
}

// ---------------------------------------------------------------------------

void GeneNetwork::load(const char *filename, GeneNetwork::format f) {
	switch (f) {
		case TSV:
			load_tsv(filename);
			break;
		case SBML:
			load_sbml(filename);
			break;
		default:
			std::cerr << "Wrong file format to load!" << std::endl;
        	exit(1);
	}
}

// ---------------------------------------------------------------------------

void GeneNetwork::sanityCheck() {
    int inSum = 0, outSum = 0;
    for (int i = 0; i < getSize(); ++i) {
        inSum += nodes_[i].getNumInputs();
        outSum += nodes_[i].getNumOutputs();
    }
    assert(inSum == outSum);
}

// ---------------------------------------------------------------------------

void GeneNetwork::load_tsv(const char *filename) {
    ifstream data_file(filename); 
    if (!data_file.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        exit(1);
    }

    std::string line, entry;
    
    // data matrix
    while (getline(data_file, line)) {
        std::istringstream ls(line);
        getline(ls, entry, '\t');
        std::string src_str = entry;
        HillGene *src = new HillGene(entry);
        //::logging::log::emit<Debug>() << entry.c_str() << "\t";
        
        getline(ls, entry, '\t');
        std::string tgt_str = entry;
		HillGene *tgt = new HillGene(entry);
        //::logging::log::emit<Debug>() << entry.c_str() << "\t";
        getline(ls, entry, '\t');
        //::logging::log::emit<Debug>() << entry.c_str() << ::logging::log::endl;
        inDegrees[tgt_str]++;
        outDegrees[src_str]++;
        degrees[tgt_str]++;
        degrees[src_str]++;
		//::logging::log::emit<Debug>() << "sign = " << entry.c_str() <<
		//::logging::log::endl;
		Edge *e;
        if (entry == tgt_str) {   // only 2 columns 
            e = new Edge(src, tgt, "+-");
        } else {
            e = new Edge(src, tgt, entry);
        }
		edges_.push_back(*e);
		if (std::find(nodes_.begin(),nodes_.end(),*src) == nodes_.end()) 
			nodes_.push_back(*src);
		if (std::find(nodes_.begin(),nodes_.end(),*tgt) == nodes_.end())
			nodes_.push_back(*tgt);
		delete src;
		delete tgt;
		delete e;
    }
    data_file.close();
    
    x_ = Vec_DP(nodes_.size());
	y_ = Vec_DP(nodes_.size());

    inDegreesVec = std::vector< std::pair<std::string,int> >(inDegrees.begin(), 
                                                             inDegrees.end());
    sort(inDegreesVec.begin(), inDegreesVec.end(), IntCmp());
        
    outDegreesVec = std::vector< std::pair<std::string,int> >(outDegrees.begin(), 
                                                              outDegrees.end());
    sort(outDegreesVec.begin(), outDegreesVec.end(), IntCmp());
    
    degreesVec = std::vector< std::pair<std::string,int> >(degrees.begin(), 
                                                           degrees.end());
    sort(degreesVec.begin(), degreesVec.end(), IntCmp());
}

// ----------------------------------------------------------------------------

/** 
 * Load a gene network from an SBML file. Overrides Structure.load(). Format must
 * be equal GeneNetwork.SBML. Note, the SBML file must be in the exact same format
 * as the SBML files produced by writeSBML(). In particular, we assume that reactions are listed
 * *ordered* as we do in writeSBML().
 * @param filename URL to the file describing the network to load
 * @param format File format (GML, DOT, etc.)
 * @throws IOException 
 */
void GeneNetwork::load_sbml(const char *filename) {
	SBMLDocument* document;
  	SBMLReader reader;

  	document = reader.readSBML(filename);

  	unsigned int errors = document->getNumErrors();
	if (errors > 0) {
        std::cerr << "Failed to open file " << filename << std::endl;
        exit(1);
    }

	Model *m = document->getModel();

	// -----------------------------------------
	// Set the network size and create the genes
	// do not count the species _void_
	int size = m->getNumSpecies() - 1;
	ListOfSpecies *species = m->getListOfSpecies();
	
	for (int g=0; g < size; g++) {
		if (species->get(g)->getId() != "_void_") {
			//HillGene hg = new HillGene(this);
			//hg.setLabel(species.get(g).getId());
			HillGene *n = new HillGene(species->get(g)->getId());
			//n.setLabel(species->get(g)->getId());
			nodes_.push_back(*n);
			delete n;
		}
	}
	
	x_ = Vec_DP(nodes_.size());
	x_ = 0;
	y_ = Vec_DP(nodes_.size());
	y_ = 0;
	
	//vector<string> parameterNames; // the names of the parameters
	//vector<double> parameterValues; // the values of the parameters
	std::map<std::string, double> params;
	std::vector<std::string> inputNodes; // the indexes of the inputs
	HillGene src, tgt;
	Parameter *param;
	
	// 2 loops for one gene: both synthesis and degradation reactions
	// (we assume that reactions are listed *ordered* as we do in writeSBML())
	//int counter = 0;
	for (unsigned int i=0; i < m->getNumReactions(); i++) {
		Reaction *re = m->getReaction(i);
		std::string id = re->getId();
		
		std::stringstream ss;
		ss << i;
		//::logging::log::emit<Debug>() << id.c_str() <<
		//		::logging::log::endl;
	
		tgt = nodes_.at(getIndexOfNode(getGeneReactantId(id)));
		//tgt->setLabel(getGeneReactantId(*re));
      	//SpeciesReference *rt = re->getReactant(0);
      	//Node *tgt = new HillGene();
      	//tgt->setLabel(rt->getSpecies());
      	//ListOfSpeciesReferences *modifiers = re->getListOfModifiers();

    	for (unsigned int j=0; j < re->getNumModifiers(); j++) {
      		ModifierSpeciesReference *md = re->getModifier(j);
      		src = nodes_.at(getIndexOfNode(md->getSpecies()));      		
      		inputNodes.push_back(src.getLabel());
      		
            // set output genes
            std::vector<std::string> outputs = src.getOutputGenes();
            outputs.push_back(tgt.getLabel());
            src.setOutputGenes(outputs);
      		
            // The edge type is unknown for now, it is initialized later
      		Edge *e = new Edge(&src, &tgt, "+-");
			edges_.push_back(*e);
			//delete src;
			delete e;
		}

      	KineticLaw *kl = re->getKineticLaw();
      		
      	for(unsigned int j=0; j < kl->getNumParameters(); j++) {
        	param = kl->getParameter(j);
			params[param->getId()] = param->getValue();
			//char buf[256];
      		//sprintf(buf, "%s\t%f", param->getId().c_str(), param->getValue());
			//::logging::log::emit<Info>() << buf <<	::logging::log::endl;
		}
		
		//::logging::log::emit<Info>() << ::logging::log::dec << params.size() <<
		//		::logging::log::endl;
		
		// in the second iteration for this gene
		if (i%2 == 1) {
			// set parameters in gene
			//tgt.initialization(params, inputNodes);
			nodes_.at(getIndexOfNode(getGeneReactantId(id))).initialization(params, inputNodes);;
			//char buf[256];
			//sprintf(buf, "%f", params["k_1"]);
			//::logging::log::emit<Info>() << buf << ::logging::log::endl;
			
			inputNodes.clear(); // don't clear because the reference was copied to the gene
			//parameterNames.clear(); // reset (they were not copied)
			//parameterValues.clear();
			params.clear();
		}
		//counter++;
	}
	//setEdgeTypesAccordingToDynamicalModel();
	//signed_ = true;
	
	//delete document;
	//delete n;
	//delete e;
}

// ----------------------------------------------------------------------------

void GeneNetwork::writeTSV(const char *filename) {
    ofstream data_file(filename); 
    if (!data_file.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        exit(1);
    }

	::logging::log::emit<Info>() << "Writing file " << filename <<
		::logging::log::endl;

    for (std::vector<Edge>::iterator it = edges_.begin(); it != edges_.end(); ++it) {
        int idx = getIndexOfNode(it->getTarget().getLabel());
        std::vector<std::string> inputGenes = nodes_[idx].getInputGenes();
        std::vector<Edge::sign> inputSigns = nodes_[idx].getInputSigns();
        std::vector<std::string>::const_iterator it_sign = 
            find(inputGenes.begin(),
                inputGenes.end(), (*it).getSource().getLabel());
        data_file << (*it).getSource().getLabel() << "\t" << 
    		(*it).getTarget().getLabel() << "\t" << 
            (*it).signToString[inputSigns[it_sign-inputGenes.begin()]] << 
            std::endl;
    }
    
    data_file.close();
}

// ----------------------------------------------------------------------------

/**
 * Save the gene network to an SBML file. If the argument is null, use the network id.
 * @param filename URL to the file describing the network to load
 * @throws IOException
 */
void GeneNetwork::writeSBML(const char *filename) {
			
	ofstream data_file(filename); 
    if (!data_file.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        exit(1);
    }
    data_file.close();
			
	::logging::log::emit<Info>() << "Writing file " << filename <<
		::logging::log::endl;
	
	SBMLDocument *sbmlDoc = new SBMLDocument(3, 1);

	Model *model = sbmlDoc->createModel();
	model->setId(id_);
	//model.getNotes ().add (comment_); // save network description
	
	int size = getSize();
	
	Compartment *comp = model->createCompartment();
  	comp->setId("cell");
	comp->setSize(1);

	std::vector<Species*> all_sp;
	Species *sp;

	for (int s=0; s < size; s++) { // save gene as species
//			species[s] = new Species(nodeIds_.get(s), nodeIds_.get(s));
		sp = model->createSpecies();
  		sp->setCompartment("cell");
  		sp->setId((nodes_.at(s)).getLabel());
  		all_sp.push_back(sp);
		//species[s].setInitialAmount(?); // maybe save the wild-type steady state?
		//model.addSpecies(species[s]);
	}
	
	// create the void species
	sp = model->createSpecies();
  	sp->setCompartment("cell");
  	sp->setId("_void_");
	sp->setInitialAmount(0);
	sp->setBoundaryCondition(true);
	sp->setConstant(true);
	all_sp.push_back(sp);
	//model.addSpecies(species[size]);


	// SET SYNTHESIS AND DEGRADATION REACTIONS FOR EVERY GENE
	for (int i=0; i<size; i++) {
		//::logging::log::emit<Info>() << ::logging::log::dec << i <<
		//::logging::log::endl;
		
		// the ID of gene i
//			String currentGeneID = nodeIds_.get(i);
		string currentGeneID = (nodes_.at(i)).getLabel();
		// The modifiers (regulators) of gene i
		std::vector<std::string> inputGenes = (nodes_.at(i)).getInputGenes();
		
		// SYNTHESIS REACTION
		std::string reactionId = currentGeneID + "_synthesis";
		Reaction *reaction = model->createReaction();
		KineticLaw *kineticLaw = reaction->createKineticLaw();
		SpeciesReference *spr;
		ModifierSpeciesReference *msr;
		reaction->setId(reactionId);
		reaction->setReversible (false);
		spr = reaction->createReactant();
  		spr->setSpecies(sp->getId());
  		spr = reaction->createProduct();
  		spr->setSpecies((all_sp.at(i))->getId());
		
		std::stringstream ss;
		ss << inputGenes.size();
		//::logging::log::emit<Debug>() << "node = " << nodes_.at(i).getLabel().c_str() << " #inputs = " << ss.str().c_str() << ::logging::log::endl;
		
		for (unsigned int r=0; r<inputGenes.size(); r++) {// set gene modifiers
//				reaction.addModifier(species[inputIndexes.get(r)]);
			//log.log(Level.INFO, "i = " + size);
			msr = reaction->createModifier();
			msr->setSpecies((all_sp.at(getIndexOfNode(inputGenes.at(r))))->getId());
		}

		//std::vector<RegulatoryModule> modules = (nodes_.at(i)).getRegulatoryModules();
		//log.log(Level.INFO, "size = " + modules.size());
		std::map<std::string, double> *params = new std::map<std::string, double>();
		(nodes_.at(i)).compileParameters(*params);
		
		//char buf[256];
		//sprintf(buf, "%f", nodes_.at(i).getDelta());
		//::logging::log::emit<Info>() << buf << ::logging::log::endl;
		//::logging::log::emit<Info>() << ::logging::log::dec << nodes_.at(i).getAlpha().size() <<
		//		::logging::log::endl;
		
		Parameter *para;
		// save gene parameters (note, the first param is the degradation rate)
		std::map<std::string, double>::iterator p = params->begin();
		//p++;
		for (; p!=params->end(); p++) {
			//if (p == params->begin()) {
			//	p++;
			//	continue;
			//}
			//::logging::log::emit<Info>() << p->first.c_str() <<
			//	::logging::log::endl;
			if (p->first != "delta") {
				para = kineticLaw->createParameter();
				para->setId(p->first);
				para->setValue(p->second);
			}
		}
		reaction->setKineticLaw(kineticLaw);
		model->addReaction(reaction);

		// DEGRADATION REACTION
		reaction = model->createReaction();
		kineticLaw = reaction->createKineticLaw();
		reactionId = currentGeneID + "_degradation";
		reaction->setId(reactionId);
		reaction->setReversible(false);
		spr = reaction->createReactant();
  		spr->setSpecies((all_sp.at(i))->getId());
  		spr = reaction->createProduct();
  		spr->setSpecies(sp->getId());

		para = kineticLaw->createParameter();
		std::map<std::string,double>::iterator it = params->find("delta");
		para->setId(it->first);
		para->setValue(it->second);
		
		reaction->setKineticLaw (kineticLaw);
		model->addReaction (reaction);
	}
	
	// PRINT FILE
	SBMLWriter sbmlWriter;
	sbmlWriter.writeSBML(sbmlDoc, filename);
	
	delete sbmlDoc;
}

// ----------------------------------------------------------------------------
	
/**
 * Return gene id of a reaction.
 * @param r Instance of Reaction
 * @return Gene id of the reaction r
 */
std::string GeneNetwork::getGeneReactantId(std::string& id) {
	
	// Example: G1_synthesis or G1_degradation -> id = "G1"
	//StringTokenizer st = new StringTokenizer(r.getId(), "_");
	//String id = st.nextToken();
	// problem: if the id also has '_', the above code fails
	//std::string id = r.getId();
	std::string toReturn = id.substr(0, id.rfind("_"));
	return toReturn;
}

// ----------------------------------------------------------------------------

int GeneNetwork::getIndexOfNode(std::string label) {

	unsigned int i;
	for (i = 0; i < nodes_.size(); ++i) {
		if ((nodes_.at(i)).getLabel() == label)
			break;
	}
	return i;
		
}

// ----------------------------------------------------------------------------

int GeneNetwork::getIndexOfNode(HillGene& gene) {

	std::vector<HillGene>::iterator it = find(nodes_.begin(), nodes_.end(), gene);
	if (it != nodes_.end())
		return it - nodes_.begin();
	else {
		::logging::log::emit<Error>() << gene.getLabel().c_str() << " not in network!" << ::logging::log::endl;
		throw Exception();
	}
		
}

// ----------------------------------------------------------------------------
	
/** Get a string with all gene labels separated by tabs (often used as header when printing data) */
std::string GeneNetwork::getHeader(bool includeProteins) {

	std::string header = "\"" + nodes_.at(0).getLabel() + "\"";
	for (unsigned int i=1; i<nodes_.size(); i++)
		header += "\t\"" + nodes_.at(i).getLabel() + "\"";
	
	if (includeProteins)
		for (unsigned int i=1; i<nodes_.size(); i++)
			header += "\t\"" + nodes_.at(i).getLabel() + "_prot\"";
	header += "\n";
	return header;
}

// ----------------------------------------------------------------------------
	
/**
 * Compute the rate of change of all state variables x (and y if translation is modelled).
 * The result is returned concatenated into one array [dxdt dydt].
 * @param xy Used to set the gene expressions x and protein concentrations y (if translation)
 * @param dxydt Variations of x and y (output)
 */
void GeneNetwork::computeDxydt(Vec_DP& xy, Vec_DP& dxydt) {
	
	bool modelTranslation = GnwSettings::Instance()->getModelTranslation();
	//double[] dxydt = new double[xy.length];
	int size = getSize();
	
	for (int i=0; i<size; i++) {
		x_[i] = xy[i];
		std::stringstream s1,s2;
		s1 << i;
		s2 << x_[i];
		//::logging::log::emit<Debug>() << "x" << s1.str().c_str() << " = " << s2.str().c_str() << ::logging::log::endl;
	}
		
	if (modelTranslation)
		for (int i=0; i<size; i++)
			y_[i] = xy[size+i];
	else
		y_ = x_;
	
	// dxydt temporarily used to store the production rates of mRNA
	computeMRnaProductionRates(dxydt);
	
	for (int i=0; i < size; i++)
		dxydt[i] = dxydt[i] - nodes_.at(i).computeMRnaDegradationRate(x_[i]);
	
	if (modelTranslation)
		for (int i=0; i<size; i++)
			dxydt[size+i] = nodes_.at(i).getMaxTranslation()*x_[i] - nodes_.at(i).computeProteinDegradationRate(y_[i]);
}

// ---------------------------------------------------------------------------

void GeneNetwork::dynamics(Vec_I_DP& xy, Vec_O_DP& dxydt, GeneNetwork& grn) {
	
	bool modelTranslation = GnwSettings::Instance()->getModelTranslation();
	//double[] dxydt = new double[xy.length];
	int size = grn.getSize();
	
	for (int i=0; i<size; i++) {
		grn.x_[i] = xy[i];
		//std::stringstream s1,s2;
		//s1 << i;
		//s2 << x_[i];
		//::logging::log::emit<Debug>() << "x" << s1.str().c_str() << " = " << s2.str().c_str() << ::logging::log::endl;
	}
	
	if (modelTranslation)
		for (int i=0; i<size; i++)
			grn.y_[i] = xy[size+i];
	else
		grn.y_ = grn.x_;
	
	// dxydt temporarily used to store the production rates of mRNA
	grn.computeMRnaProductionRates(dxydt);
	
	for (int i=0; i < size; i++)
		dxydt[i] = dxydt[i] - grn.nodes_.at(i).computeMRnaDegradationRate(grn.x_[i]);
	
	if (modelTranslation)
		for (int i=0; i<size; i++)
			dxydt[size+i] = grn.nodes_.at(i).getMaxTranslation()*grn.x_[i] - 
                grn.nodes_.at(i).computeProteinDegradationRate(grn.y_[i]);
}

// ---------------------------------------------------------------------------

int GeneNetwork::gene_network_dynamics(double, const double xy[], double dxydt[], void *param) { 
	
	GeneNetwork *grn = (GeneNetwork*) param;
  	bool modelTranslation = GnwSettings::Instance()->getModelTranslation();
	//double[] dxydt = new double[xy.length];
	int size = grn->getSize();
	
	for (int i=0; i<size; i++) {
        if (xy[i] >= 0)
		    grn->x_[i] = xy[i];
        else
            grn->x_[i] = 0;
		std::stringstream s1,s2;
		s1 << i;
		s2 << grn->x_[i];
		//::logging::log::emit<Debug>() << "x" << s1.str().c_str() << " = " << s2.str().c_str() << ::logging::log::endl;
	}
	
	if (modelTranslation)
		for (int i=0; i<size; i++)
			grn->y_[i] = xy[size+i];
	else
		grn->y_ = grn->x_;
	
	// dxydt temporarily used to store the production rates of mRNA
	Vec_DP dxydt_v(size); 
	grn->computeMRnaProductionRates(dxydt_v);
	
	for (int i=0; i < size; i++)
		dxydt_v[i] = dxydt_v[i] - grn->getNodes().at(i).computeMRnaDegradationRate(grn->x_[i]);
		
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
// PRIVATE FUNCTIONS
	
/**
 * Compute the production rate for all the genes. This function is identical to
 * the one below, except that the argument is passed as double[] instead of
 * DoubleMatrix1d. We need both because the SDE implementation works with
 * DoubleMatrix1d and the ODE implementation works with arrays.
 * NOTE: x_ must be set before calling this function.
 * @param productionRates Return the production rates of all the genes
 */
inline
void GeneNetwork::computeMRnaProductionRates(Vec_DP& productionRates) {
	
	int size = getSize();
	
	for (int i=0; i < size; i++) {
		std::vector<std::string> inputGenes = nodes_.at(i).getInputGenes();
		Vec_INT inputIndex(inputGenes.size());
		
		std::stringstream ss, ssx;
		//ss << nodes_.at(i).getRegulatoryModules().size();
        ssx << x_[i];
        if (x_[i] < 0 || std::isinf(x_[i]))
            ::logging::log::emit<Debug>() << "node = " << 
                nodes_.at(i).getLabel().c_str() << " x = " << 
                ssx.str().c_str() << ::logging::log::endl;
		
		if (GnwSettings::Instance()->getModelTranslation()) {
			for (int j = 0; j < inputIndex.size(); ++j)
				inputIndex[j] = getIndexOfNode(inputGenes.at(j));
			productionRates[i] = nodes_.at(i).computeMRnaProductionRate(inputIndex, y_);
		}
		else {
			for (int j = 0; j < inputIndex.size(); ++j)
				inputIndex[j] = getIndexOfNode(inputGenes.at(j));
			productionRates[i] = nodes_.at(i).computeMRnaProductionRate(inputIndex, x_);
		}
	}
}

// ---------------------------------------------------------------------------

/**
 * Compute the production rate for all the genes. This function is identical to
 * the one above, except that the argument is passed as DoubleMatrix1d instead of
 * double[]. We need both because the SDE implementation works with
 * DoubleMatrix1d and the ODE implementation works with arrays.
 * NOTE: x_ must be set before calling this function.
 * @param productionRates Return the production rates of all the genes
 */
void GeneNetwork::computeInputs(Vec_I_DP& x, Vec_O_DP& productionRates) {
	
	int size = getSize();
	
	for (int i=0; i < size; i++) {
		//if (GnwSettings.getInstance().getModelTranslation())
		//	productionRates.set(i, ((Gene)nodes_.get(i)).computeMRnaProductionRate(i, y_));
		//else
        std::stringstream ss;
        ss << x[i];
		if (x[i] < 0)
            ::logging::log::emit<Debug>() << "node = " << nodes_.at(i).getLabel().c_str() << " x = " << ss.str().c_str() << ::logging::log::endl;
		
        std::vector<std::string> inputGenes = nodes_.at(i).getInputGenes();
		Vec_INT inputIndex(inputGenes.size());
		for (int j = 0; j < inputIndex.size(); ++j)
			inputIndex[j] = getIndexOfNode(inputGenes.at(j));
		productionRates[i] = nodes_.at(i).computeMRnaProductionRate(inputIndex, x);
		
		std::stringstream ss1,ss2,ss3;
		ss1 << i;
		ss2 << productionRates[i];
		ss3 << x[i];
		//::logging::log::emit<Debug>() << ss1.str().c_str() << " " <<  ss2.str().c_str() << " " << ss3.str().c_str() << ::logging::log::endl;
	}
}

// ---------------------------------------------------------------------------

void GeneNetwork::setRank(const char *filename) {
	ifstream data_file(filename); 
    if (!data_file.is_open()) {
        ::logging::log::emit<Error>() << "Failed to open file " << filename << ::logging::log::endl;
		throw Exception();
    }

    std::string line, entry;
    
    // ranked list
    int counter = 1;
    while (getline(data_file, line)) {
        std::istringstream ls(line);
        getline(ls, entry, '\n');
		rank_.push_back(nodes_.at(getIndexOfNode(entry)));
		counter++;
    }
    data_file.close();
}

// ---------------------------------------------------------------------------

std::vector<HillGene>& GeneNetwork::getRank() {
	return rank_;
}

// ----------------------------------------------------------------------------

/**
 * Return a list of all genes that perturbe the gene given in parameter.
 * @param gene Gene perturbed
 * @return Instance of all genes that perturbed the gene given in parameter.
 */
void GeneNetwork::getInputs(HillGene& gene, std::vector<std::string>& inputs) {

	for (unsigned int i=0; i < edges_.size(); i++) {
		if (edges_.at(i).getTarget().getLabel() == gene.getLabel()) {
			std::string source = edges_.at(i).getSource().getLabel();
			inputs.push_back(nodes_.at(getIndexOfNode(source)).getLabel());
		}
	}
	
}

// ----------------------------------------------------------------------------

/**
 * Return a list of genes that are output target of the given gene.
 * @param gene The gene of interest
 * @param inputs Return value, list of genes that are output of the gene of
 * interest
 */
void GeneNetwork::getOutputs(HillGene& gene, std::vector<std::string>& outputs) {

	for (unsigned int i=0; i < edges_.size(); i++) {
		if (edges_.at(i).getSource().getLabel() == gene.getLabel()) {
			std::string target = edges_.at(i).getTarget().getLabel();
			outputs.push_back(target);
		}
	}
	
}

// ----------------------------------------------------------------------------

void GeneNetwork::getInputSigns(HillGene& gene, std::vector<Edge::sign>& inputSigns) {

	for (unsigned int i=0; i < edges_.size(); i++) {
		if (edges_.at(i).getTarget().getLabel() == gene.getLabel()) {
			//std::string source = edges_.at(i).getSource().getLabel();
			inputSigns.push_back(edges_.at(i).getSign());
            //::logging::log::emit<Debug>() << 
            //    edges_[i].signToString[edges_[i].getSign()].c_str() << 
            //    ::logging::log::endl;
		}
	}
	
}

// ----------------------------------------------------------------------------
	
/**
 * Based on W_, determine for each gene which are its inputs
 */
void GeneNetwork::initializeInputWiring() {
	
	int size = getSize();
	
	for (int i=0; i < size; i++) {
		std::vector<std::string> inputs; 
		getInputs(nodes_.at(i), inputs);
		nodes_.at(i).setInputGenes(inputs);
		
		std::vector<std::string> outputs; 
		getOutputs(nodes_.at(i), outputs);
		nodes_.at(i).setOutputGenes(outputs);
		
		std::vector<Edge::sign> inputSigns; 
		getInputSigns(nodes_.at(i), inputSigns);
		nodes_.at(i).setInputSigns(inputSigns);
		
	}
}

// ----------------------------------------------------------------------------
	
/**
 * Random initialization of the dynamical model (based on the fixed topology defined by W_)
 */
void GeneNetwork::randomInitialization() {
	
	::logging::log::emit<Info>() << "Random initialization ..." <<  ::logging::log::endl;

	// set the inputs for the genes based on W_
	initializeInputWiring();
	int size = nodes_.size();

	//HillGene gene;
	for (int i=0; i < size; i++) {
		nodes_.at(i).randomKineticInitialization();
	}
    std::stringstream ss1, ss2;
	for (int i=0; i < size; i++) {
        ss1 << i + 1;
        ss2 << nodes_[i].getNumInputs();
		//::logging::log::emit<Debug>() << ss1.str().c_str() << "\t" << 
        //    nodes_.at(i).getLabel().c_str() << "\t" << ss2.str().c_str()
        //    << ::logging::log::endl;
		nodes_.at(i).randomStructureInitialization();
		/*
		std::stringstream ss;
		if (nodes_.at(i).getRegulatoryModules().size() > 0) {
			ss << nodes_.at(i).getLabel() << " " << nodes_.at(i).getRegulatoryModules().at(0).getNumInputs();
			::logging::log::emit<Debug>() << ss.str().c_str() << ::logging::log::endl;
		}
		*/
		//log.log(Level.INFO, gene.getLabel() + ": " + gene.toString());
	}
	// If the network was unsigned before the initialization, we have to set the types
	// of the edges (enhancing / inhibiting) accordingly. If the network was signed,
	// this just sets the same values that are already there.
	//setEdgeTypesAccordingToDynamicalModel();
	//signed_ = true;
}

// ----------------------------------------------------------------------------

/**
 * Set the type (enhancing, inhibiting, dual, or unknown) of all edges according
 * to the dynamically model.
 */
/*
void GeneNetwork::setEdgeTypesAccordingToDynamicalModel() {
	
	for (unsigned int i=0; i<nodes_.size(); i++) {
		nodes_.at(i).setInputEdgeTypesAccordingToDynamicalModel();
	}
}
*/

// ---------------------------------------------------------------------------

void GeneNetwork::getEdge(HillGene& src, HillGene& tgt, Edge& e) {
	for (unsigned int i = 0; i < edges_.size(); ++i) {
		if (edges_.at(i).getSource() == src && edges_.at(i).getTarget() == tgt) {
			e = edges_.at(i);
			return;
		}
	}
	::logging::log::emit<Info>() << "Edge " << src.getLabel().c_str() << "->" << tgt.getLabel().c_str() << " not found!" <<  ::logging::log::endl;
}

// ============================================================================
	
/**
 * implement the pruning
 */
void GeneNetwork::prune(std::vector<HillGene>& genes) {
	
	//log.log(Level.INFO, "edge_size = " + edges_.size());
	std::vector<Edge> edge_clone;
	for (unsigned int iNode = 0; iNode < genes.size(); iNode++) {
		::logging::log::emit<Info>() << "Pruning " << genes.at(iNode).getLabel().c_str() <<  ::logging::log::endl;
	
		//edge_clone = new ArrayList<Edge>();
		for (unsigned int iEdge = 0; iEdge < edges_.size(); iEdge++) {
			//log.log(Level.INFO, "edge " + Integer.toString(iEdge) + " " + edges_.get(iEdge).toString() + " " + genes.get(iNode));
			if (edges_.at(iEdge).getSource() == genes.at(iNode) || 
				edges_.at(iEdge).getTarget() == genes.at(iNode)) {
				//edges_.erase(edges_.begin()+iEdge);
		        ::logging::log::emit<Info>() << "\tPruning " << edges_.at(iEdge).getSource().getLabel().c_str() << " -> " << edges_.at(iEdge).getTarget().getLabel().c_str() <<  ::logging::log::endl;
				nodes_.at(getIndexOfNode(edges_.at(iEdge).getTarget().getLabel())).pruneInput(nodes_.at(getIndexOfNode(edges_.at(iEdge).getSource().getLabel())));
                //std::stringstream ss;
                //ss << nodes_.at(getIndexOfNode(edges_.at(iEdge).getTarget().getLabel())).getRegulatoryModules().size();
		        //::logging::log::emit<Debug>() << nodes_.at(getIndexOfNode(edges_.at(iEdge).getTarget().getLabel())).getLabel().c_str() << " input size " << ss.str().c_str() << ::logging::log::endl;
			} else {
				edge_clone.push_back(edges_.at(iEdge));
			}
		}
		edges_ = edge_clone;
	}
	edges_ = edge_clone;
	//log.log(Level.INFO, "edge_size = " + edges_.size());
	
	std::vector<HillGene> node_clone;
	for (unsigned int iEdge = 0; iEdge < edges_.size(); iEdge++) {
		if (find(node_clone.begin(),node_clone.end(),edges_.at(iEdge).getSource()) == node_clone.end())
			node_clone.push_back(nodes_.at(getIndexOfNode(edges_.at(iEdge).getSource().getLabel())));
		if (find(node_clone.begin(),node_clone.end(),edges_.at(iEdge).getTarget()) == node_clone.end())
			node_clone.push_back(nodes_.at(getIndexOfNode(edges_.at(iEdge).getTarget().getLabel())));
	}
	nodes_ = node_clone;
	//log.log(Level.INFO, "size = " + nodes_.size());
    //std::stringstream ss;
    //ss << nodes_.at(getIndexOfNode("G37")).getRegulatoryModules().size();
    //::logging::log::emit<Debug>() << "G37 module size " << ss.str().c_str() << ::logging::log::endl;
	
	//setSize(getSize());
	
	initializeInputWiring();
	
	// remove k,n,alpha
	std::vector<RegulatoryModule> modules;
	std::vector<HillGene> inputs;
	//ArrayList<String> label;
	HillGene src;
	Vec_DP k, n;
	for (unsigned int iNode = 0; iNode < nodes_.size(); iNode++) {
		/*
		modules = nodes_.at(iNode).getRegulatoryModules();
		for (int iModule = 0; iModule < modules.size(); iModule++) {
			k = modules.at(iModule).getK();
			n = modules.at(iModule).getN();
			inputs = modules.at(iModule).getInputs(); 
			for (int iInput = 0; iInput < inputs.size(); iInput++) {
				label = new ArrayList<String>();
				label.add(inputs.get(iInput).getLabel());
				src = nodes_.get(getNodeIndexesFromNodeLabels(label).get(0));
				if (!containsEdge(getEdge(src, nodes_.get(iNode)))) {
					n[iInput] = 0;
				}
			}
			modules.get(iModule).setK(k);
		}
		*/
		//nodes_.at(iNode).randomInitializationOfAlpha();
	}
	
	//return *this;
		
}

