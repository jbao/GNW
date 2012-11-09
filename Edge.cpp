#include "Edge.h"

Edge::Edge(const Node *src, const Node *tgt, std::string s) : src_(*src), tgt_(*tgt) {
	//stringToSign = new std::map<std::string, sign>();
	//signToString = new std::map<sign, std::string>();
	stringToSign["+"] = ACTIVATOR;
	stringToSign["-"] = INHIBITOR;
	stringToSign["+-"] = UNKNOWN;
	signToString[ACTIVATOR] = "+";
	signToString[INHIBITOR] = "-";
	signToString[UNKNOWN] = "+-";
	
	sign_ = stringToSign[s];
}

Edge::~Edge() {
	//delete stringToSign;
	//delete signToString;
}