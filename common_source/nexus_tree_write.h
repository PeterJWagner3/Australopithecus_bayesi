#ifndef NEXUS_TREE_WRITE_H
#define NEXUS_TREE_WRITE_H

#include <stdio.h>
#include "node.h"
#include "tree_fns.h"

class WriteNexus
	{
	public:
		void WriteTreeBlock(NODE *tree, int *ntaxa);
		void WriteNexusTreeFile (NODE **, int *, int *, char **, char *);
		void WriteNexusTreeFile (NODE **, int *, int *);
		void WriteNexusTreeFile (NODE *, int *);
		void WriteNexusTreeWithBL (FILE *nexus, NODE *node);
		void WriteNexusTree (FILE *, NODE *);

	private:
		FILE *nexus;

	};

#endif

