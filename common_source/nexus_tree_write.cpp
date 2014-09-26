#include "nexus_tree_write.h"

void WriteNexus::WriteTreeBlock(NODE *tree, int *ntaxa)
{
	int i;
	NODE *root;

	fprintf (nexus, "BEGIN TREES;\n\n");
	fprintf (nexus, "\tTRANSLATE\n");
	for (i=1;i<*ntaxa;++i) fprintf (nexus, "\t\t%d\tTaxon%d,\n", i, i);
	fprintf (nexus, "\t\t%d\tTaxon%d\n;\n\n", *ntaxa, *ntaxa);
	root=GetRoot(tree, ntaxa);
	fprintf (nexus, "\tTREE\tTree_%d = \t[&R] ", i);
	WriteNexusTree(nexus, root);
	fprintf (nexus, ";\n");
	fprintf (nexus, ";\n\nEND;\n\n");
}

void WriteNexus::WriteNexusTreeFile(NODE *tree, int *ntaxa)
{
	nexus=fopen("treefile.nex", "w");
	fprintf (nexus, "#NEXUS\n\n");
	WriteTreeBlock(tree, ntaxa);			
	fclose (nexus);
}

void WriteNexus::WriteNexusTreeFile (NODE **tree_list, int *ntrees, int *ntaxa)
{
	FILE *nexus;
	NODE *root;
	int i;

	nexus=fopen("treefile.nex", "w");
	fprintf (nexus, "#NEXUS\n\n");
	fprintf (nexus, "BEGIN TREES;\n\n");
	fprintf (nexus, "\tTRANSLATE\n");
	for (i=1;i<*ntaxa;++i) fprintf (nexus, "\t\t%d\tTaxon%d,\n", i, i);
	fprintf (nexus, "\t\t%d\tTaxon%d\n;\n\n", *ntaxa, *ntaxa);
	for (i=1;i<=*ntrees;++i)
		{
		root=GetRoot(tree_list[i], ntaxa);
		fprintf (nexus, "\tTREE\tTree_%d = \t[&R]\t", i);
		WriteNexusTree(nexus, root);
		fprintf (nexus, ";\n");
		}
	fprintf (nexus, ";\n\nEND;\n\n");
	fclose (nexus);
}

void WriteNexus::WriteNexusTreeFile (NODE **tree_list, int *ntrees, int *ntaxa, char **taxon_label, char *filename)
{
	FILE *nexus;
	NODE *root;
	int i;

	nexus=fopen(filename, "w");
	fprintf (nexus, "#NEXUS\n\n");
	fprintf (nexus, "BEGIN TREES;\n\n");
	fprintf (nexus, "\tTRANSLATE\n");
	for (i=1;i<*ntaxa;++i) fprintf (nexus, "\t\t%d\t%s,\n", i, taxon_label[i]);
	fprintf (nexus, "\t\t%d\t%s\n;\n\n", *ntaxa, taxon_label[*ntaxa]);
	for (i=1;i<=*ntrees;++i)
		{
		root=GetRoot(tree_list[i], ntaxa);
		fprintf (nexus, "\tTREE\tTree_%d = \t[&R]\t", i);
		WriteNexusTree(nexus, root);
		fprintf (nexus, ";\n");
		}
	fprintf (nexus, ";\n\nEND;\n\n");
	fclose (nexus);
}

void WriteNexus::WriteNexusTree (FILE *nexus, NODE *node) {
	int i;
	NODE *desc;
	NODE *sis;

	desc=node->left;
	sis=node->right;
	
	if (!desc->isAnc && !sis->isAnc) {
		for (i=0; i<=1; ++i) {
			if (i) {
				desc=node->right;
				sis=node->left;
			}
		if (!i) fprintf(nexus, "(");
		if (desc->tip) fprintf(nexus, "%d", desc->label);
		else WriteNexusTree(nexus, desc);
//			if (desc->bl) fprintf (nexus, ":%f", *desc->bl);
		if (!i) fprintf(nexus, ",");
		else fprintf(nexus, ")");
		}
	}
	else {
		if (desc->isAnc) {
			fprintf(nexus, "(");
			if (sis->tip) fprintf(nexus, "%d", sis->label);
			else WriteNexusTree(nexus, sis);
//			if (desc->bl) fprintf (nexus, ":%f", *desc->bl);
			fprintf(nexus, ")%d", desc->label);
		}
		else {
			fprintf(nexus, "(");
			if (desc->tip) fprintf(nexus, "%d", desc->label);
			else WriteNexusTree(nexus, desc);
//			if (desc->bl) fprintf (nexus, ":%f", *desc->bl);
			fprintf(nexus, ")%d", sis->label);
		}
	}
}

void WriteNexus::WriteNexusTreeWithBL (FILE *nexus, NODE *node) {
	// it appears that it is impossible in NEXUS format to have both assigned ancestors AND
	// branch lengths (e.g., "1(2):1.4"

	int i;
	NODE *desc;
	NODE *sis;

	desc=node->left;
	sis=node->right;
	
//	if (!desc->isAnc && !sis->isAnc) {
		for (i=0; i<=1; ++i) {
			if (i) {
				desc=node->right;
				sis=node->left;
			}
		if (!i) fprintf(nexus, "(");
		if (desc->tip) 
			{
			fprintf(nexus, "%d", desc->label);
			fprintf (nexus, ":%f", *desc->bl);
			}
		else WriteNexusTreeWithBL(nexus, desc);
		if (!i) fprintf(nexus, ",");
		else 
			{
			fprintf(nexus, ")");
			fprintf (nexus, ":%f", *node->bl);
			}
		}
/*	}
	else {
		if (desc->isAnc) {
			fprintf(nexus, "(");
			if (sis->tip) {
				fprintf(nexus, "%d", sis->label);
				fprintf (nexus, ":%f", *sis->bl);
			}
		else WriteNexusTreeWithBL(nexus, sis);
		fprintf(nexus, ")%d", desc->label);
		fprintf (nexus, ":%f", *desc->bl);
		}
		else {
			fprintf(nexus, "(");
			if (desc->tip) {
				fprintf(nexus, "%d", desc->label);
				fprintf (nexus, ":%f", *desc->bl);
			}
		else WriteNexusTreeWithBL(nexus, desc);
		fprintf(nexus, ")%d", sis->label);
		fprintf (nexus, ":%f", *sis->bl);
		}
	}
*/
}

