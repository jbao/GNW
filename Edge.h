#ifndef EDGE_H
#define EDGE_H

#include <map>
#include "Node.h"

class Edge {

public:
	enum sign {ACTIVATOR, INHIBITOR, UNKNOWN};
	std::map<std::string, sign> stringToSign;
	std::map<sign, std::string> signToString;
	Edge(const Node *src, const Node *tgt, std::string s);
	~Edge();
	inline bool operator== (Edge &rhs) const {
        return src_ == rhs.getSource() && tgt_ == rhs.getTarget();
    }
	
	Node getSource() { return src_; }
	void setSource(Node src) { src_ = src; }
	
	Node getTarget() { return tgt_; }
	void setTarget(Node tgt) { tgt_ = tgt; }
	
	sign getSign() { return sign_; }
	void setSign(sign s) { sign_ = s; }
	
private:
	Node src_, tgt_;
	sign sign_;

};

#endif
