#include "node.h"

using namespace std;

NODE::NODE() {
	left=0;
	right=0;
	anc=0;
	sister=0;
	tip=false;
	isAnc=false;
	label=0;
	bl=0;
	richness=0;
}

NODE::~NODE() {}

void NODE::Fix ()
	{
	if (tip && !sister) isAnc=true;
	else if (tip && !sister->isAnc) isAnc=true;
	else if (tip && sister->isAnc) cout << "\n***\tAttempted to fix a tip ("<< label << ") that was already fixed (" << sister->label << ")\n";
	else cout << "\n***\tAttempted to fix an internal node " << label << "\n";
	}

void NODE::Unfix () { isAnc=false; }

