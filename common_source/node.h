#ifndef NODE_H
#define NODE_H

#include <iostream>

class NODE {
	public:
		NODE();
		~NODE();
		NODE *left;			// left descendant node
		NODE *right;		// right descendant node
		NODE *anc;			// ancestral node
		NODE *sister;		// sister node
		bool tip;			// if node is a tip (OTU)
		bool isAnc;			// if node is fixed as an ancestor
		int label;			// integer label of the node
		float *bl;			// branch length below current node
		int *richness;		// number of OTUs in the clade

		void Fix();
		void Unfix();
		
//	private:
};

#endif

