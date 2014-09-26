#include "data_field.h"

using namespace std;

DATA::DATA() {
	ntaxa = new int(0);
	nchrs = new int(0);
	DNA=0;
	stratchr=0;
	matrix=0;
	nstates=0;
	chrtype=0;
	polymorph=0;
	hasPolymorph=false;

	fint=0;
	lint=0;
	hasStratChr=false;
	FAChrPresent=false;
	LAChrPresent=false;
	MaxInt=0;
}

DATA::DATA(int no_taxa, int no_chrs) {
	ntaxa = new int(no_taxa);
	nchrs = new int(no_chrs);
	DNA=0;
	stratchr=0;
	matrix=0;
	nstates=0;
	chrtype=0;
	polymorph=0;
	hasPolymorph=false;

	fint=0;
	lint=0;
	hasStratChr=false;
	FAChrPresent=false;
	LAChrPresent=false;
	MaxInt=0;

	AllocateMatrix();
//	AllocateChrType();
}

//DATA::~DATA() {}
DATA::~DATA() {
	if (matrix) free_matrix(matrix, 2**ntaxa);
	if (polymorph) free_matrix(polymorph, 2**ntaxa);
	if (ntaxa) delete ntaxa;
	if (nchrs) delete nchrs;
	if (nstates) delete [] nstates;
	if (chrtype) delete [] chrtype;
	if (fint) delete [] fint;
	if (lint) delete [] lint;
}

void DATA::AllocateMatrix() {
	if (!nchrs) cout << "\n\n***\tCan't allocate matrix without nchrs\n";
	else if (!ntaxa) cout << "\n\n***\tCan't allocate matrix without ntaxa\n";
	else {
		matrix=ullmatrix(2**ntaxa, *nchrs+1);
		polymorph=bmatrix(2**ntaxa, *nchrs+1);
	}
	AllocateTaxonLabels();
	AllocateChrType();
	AllocateChrWeights();
	AllocateIncludedChrs();
}

void DATA::DefaultTaxonLabels() {
	int i;
 	ostringstream ss;
	
	for (i=1;i<=*ntaxa;++i) {
		ss<<"Taxon_"<<i<<"\0";
		strcpy(taxon_label[i], ss.str().c_str());
		ss.str("");
	}
}

void DATA::AddChrsToMatrix(int NumNewChrs) {
	int i, j;
	unsigned long long **old_matrix;
	unsigned long long **new_matrix;
	bool **old_polymorph;
	bool **new_polymorph;
	
	old_matrix=matrix;
	new_matrix=ullmatrix(2**ntaxa, *nchrs+NumNewChrs+1);
	old_polymorph=polymorph;
	new_polymorph=bmatrix(2**ntaxa, *nchrs+NumNewChrs+1);

	for (i=0;i<2**ntaxa;++i) {
		for (j=0;j<=*nchrs;++j) {
			new_matrix[i][j]=matrix[i][j];
			new_polymorph[i][j]=polymorph[i][j];
		}
	}

	matrix=new_matrix;
	polymorph=new_polymorph;
	free_matrix(old_matrix, 2**ntaxa);
	free_matrix(old_polymorph, 2**ntaxa);
	
	*nchrs+=NumNewChrs;

	if (chrtype) delete [] chrtype;
	if (weight) delete [] weight;
	if (included) delete [] included;

	if (nstates) delete [] nstates;
	if (infchr) delete [] infchr;	
	
	AllocateChrType();
	AllocateChrWeights();
	AllocateIncludedChrs();
	
	// Should be performed after new character is filled
	MakeNStates();
//	MakeInfChrs();
}

void DATA::AllocateTaxonLabels() {
	if (!ntaxa) cout << "\n\n***\tCan't allocate taxon labels without ntaxa\n";
	else taxon_label=smatrix(*ntaxa+1, MAXLABELSIZE);
}

void DATA::AllocateChrType() {
	if (!nchrs) cout << "\n\n***\tCan't allocate type array without nchrs\n";
	else {
		chrtype= new int [*nchrs+1];
		for (int i=0;i<=*nchrs;++i) chrtype[i]=0;
	}
}

void DATA::AllocateChrWeights() {
	if (!nchrs) cout << "\n\n***\tCan't allocate type array without nchrs\n";
	else {
		weight = new float [*nchrs+1];
		for (int i=0;i<=*nchrs;++i) weight[i]=1;
	}
}

void DATA::AllocateIncludedChrs() {	
	if (!nchrs) cout << "\n\n***\tCan't allocate type array without nchrs\n";
	else {
		included=new bool[*nchrs+1];
		for (int i=1;i<=*nchrs;++i) included[i]=true;
	}
}
	
void DATA::SetDNA (bool i) { DNA=i; }
void DATA::SetNTaxa (int i) { *ntaxa=i; }
void DATA::SetNChrs (int i) { *nchrs=i; }
void DATA::SetStratChr(int i) { stratchr=i; }
void DATA::SetMatrix (int i, int j, unsigned long long value) { matrix[i][j]=value; }
void DATA::SetNStates (int requested_chr, int i) { nstates[requested_chr]=i; }
void DATA::SetChrType (int requested_chr, int chr_type) { chrtype[requested_chr]=chr_type; }
void DATA::SetFInt (int taxon, int interval) { fint[taxon]=interval; }
void DATA::SetLInt (int taxon, int interval) { lint[taxon]=interval; }
void DATA::SetFInt (NODE *node, int interval) { fint[node->label]=interval; }
void DATA::SetLInt (NODE *node, int interval) { lint[node->label]=interval; }
void DATA::ExcludeChr(int chr_number) { included[chr_number]=false; }
void DATA::IncludeChr(int chr_number) { included[chr_number]=true; }
void DATA::SetInfChr (int chr_number, bool i) { infchr[chr_number]=i; }
void DATA::SetChrWeight (int chr_number, float i) { weight[chr_number]=i; }
void DATA::SetPolymorph (int taxon, int chr_number, bool setting) { 
	polymorph[taxon][chr_number]=setting; 
	// the following code sets resets the polymorphic indicator (polymorph[0][chr]) and hasPolymorph as the polymorphic chr thing is changed
	
	if (chrtype[chr_number]==4) {
		for (int i=1;i<2**ntaxa;++i) polymorph[i][chr_number]=false;
	}
	if (taxon<=*ntaxa) {
		if (!setting) {
			int i=1;
			polymorph[0][chr_number]=false;
			while (!polymorph[0][chr_number] && i<=*ntaxa) {
				if (polymorph[i][chr_number]) polymorph[0][chr_number]=true;
				++i;
			}
			if (!polymorph[0][chr_number]) { // if we've turned off polymorph on this chr check the other chrs.
				i=1;
				hasPolymorph=false;
				while (!hasPolymorph && chrtype[chr_number]!=4 && i<=*nchrs) {
					if (polymorph[0][chr_number]) hasPolymorph=true;
					++i;
					}
				}
			}
		else if (chrtype[chr_number]!=4) {	
			polymorph[0][chr_number]=true;
			hasPolymorph=true;
		}
	}
}

unsigned long long **DATA::PassMatrix() const { return matrix; }
unsigned long long **DATA::PassBootMatrix() const { return boot_matrix; }
bool DATA::IsDNA() const { return DNA; }
int *DATA::GetNTaxa ()  const{ return ntaxa; }
int *DATA::GetNChrs ()  const{ return nchrs; }
int DATA::GetStratChr ()  const{ return stratchr; }
unsigned long long *DATA::GetMatrix (int i, int j)  const{ return &matrix[i][j]; }
int *DATA::GetNStates (int requested_chr) const { return &nstates[requested_chr]; }
int *DATA::GetMaxObsState () { return &maxobsstate; }
int *DATA::GetChrType (int requested_chr) const { return &chrtype[requested_chr]; }
int DATA::GetMaxInt() const { return MaxInt; }
int *DATA::GetFInt (int taxon) const { return &fint[taxon]; }
int *DATA::GetLInt (int taxon) const { return &lint[taxon]; }
int *DATA::GetFInt (NODE *node) const { return &fint[node->label]; }
int *DATA::GetLInt (NODE *node) const { return &lint[node->label]; }
int *DATA::GetPattern(int chr_number) const { return &pattern[chr_number]; }
bool DATA::GetInc(int chr_number) const { return included[chr_number]; }
bool DATA::GetInfChr(int chr_number) const { return infchr[chr_number]; }
float *DATA::GetChrWeight(int chr_number) const { return &weight[chr_number]; }
bool DATA::GetPolymorph(int taxon, int chr_number) const {
	if (!polymorph) {
		cout << "\n\n***\tPolymorphic matrix has not yet beend constructed\n";
		return 0;
	}
	else return polymorph[taxon][chr_number]; 
}

void DATA::MakeMatrixFromFile(char *filename) {
	int i, j;
	FILE *nexus_file;

	AllocateMatrix();
	if ((nexus_file = fopen(filename, "r")) == NULL) {
		cout <<  "Error opening data file.\n";
		exit (0);
	}

	for (i=1;i<=*ntaxa;++i) {
		GetTaxonNameFromFile(nexus_file, &i);
		j=1;
		while (GetStateFromFile(nexus_file, &i, &j)) ++j;
	}
		
	fclose (nexus_file);
	MakeNStates();
//	MakeInfChrs();
}

void DATA::MakeMatrixFromFile(char *filename, long position) {
	int i, j;
	FILE *nexus_file;

	AllocateMatrix();
	if ((nexus_file = fopen(filename, "r")) == NULL) {
		cout <<  "Error opening data file.\n";
		exit(0);
	}

	fseek (nexus_file,position,SEEK_SET);
	for (i=1;i<=*ntaxa;++i) {
		GetTaxonNameFromFile(nexus_file, &i);
		j=1;
		while (GetStateFromFile(nexus_file, &i, &j)) ++j;
	}
		
	fclose (nexus_file);
	MakeNStates();
//	MakeInfChrs();
}


void DATA::GetTaxonNameFromFile(FILE *nexus_file, int * const taxon) {
	char inchar;
	int i;

	inchar=fgetc(nexus_file);
	while (inchar==' ' || inchar=='\t') inchar=fgetc(nexus_file);
	i=0;
	if (inchar==39) {			// ASCII 39 is the single quote
		taxon_label[*taxon][i]=inchar;
		inchar=fgetc(nexus_file);
		while (inchar!=39) {
			++i;
			taxon_label[*taxon][i]=inchar;
			inchar=fgetc(nexus_file);
		} 
	++i;
	taxon_label[*taxon][i]=inchar;
	++i;
	}
	else while (inchar!=' ' && inchar!='\t') {
		taxon_label[*taxon][i]=inchar;
		inchar=fgetc(nexus_file);
		++i;
	}
	taxon_label[*taxon][i]='\0';
	return;
}

bool DATA::GetStateFromFile(FILE *nexus_file, int *taxon, int *chr_number) {
	char inchar;
	int obs_state;
	
	inchar=fgetc(nexus_file);
	while (inchar==' ' || inchar=='\t') inchar=fgetc(nexus_file);
	if (inchar=='\n' || inchar==EOF || inchar=='\r' || inchar==0x0d || inchar==0x0a) return false;
	else if (inchar!='?' && inchar!='-') {
		if (inchar=='(') {
			hasPolymorph=true;
			polymorph[*taxon][*chr_number]=true;
			polymorph[0][*chr_number]=true;
			inchar=fgetc(nexus_file);
			if (inchar==',') inchar=fgetc(nexus_file);
			while (inchar!=')') {
				if (!ChrToState(&inchar, &obs_state)) {
					cout << "\n\n\t\t****\tError converting state "<<inchar<<" to integer in chr "<< *chr_number << " of Taxon " << *taxon << endl;
					exit(0);
				}
				else matrix[*taxon][*chr_number]+=pow(2, obs_state);
				inchar=fgetc(nexus_file);
				if (inchar==',') inchar=fgetc(nexus_file);
			}
		}
		else if (inchar=='{') {
			polymorph[*taxon][*chr_number]=false;
			inchar=fgetc(nexus_file);
			if (inchar==',') inchar=fgetc(nexus_file);
			while (inchar!='}') {
				if (!ChrToState(&inchar, &obs_state)) {
					cout << "\n\n\t\t****\tError converting state "<<inchar<<" to integer in chr "<< *chr_number << " of Taxon " << *taxon << endl;
					exit(0);
				}
				else matrix[*taxon][*chr_number]+=pow(2,obs_state);
				inchar=fgetc(nexus_file);
				if (inchar==',') inchar=fgetc(nexus_file);
			}
		}
		else {
			polymorph[*taxon][*chr_number]=false;
			if (!ChrToState(&inchar, &obs_state)) {
				cout << "\n\n\t\t****\tError converting state "<<inchar<<" to integer in chr "<< *chr_number << " of Taxon " << *taxon << endl;
				exit(0);
			}
			else matrix[*taxon][*chr_number]+=pow(2,obs_state);
		}
	}
//	delete inchar;
	return true;

}

bool DATA::ChrToState(char *inchar, int *obs_state) {
	stringstream ss;
//	int gobots;
	
	ss.str(inchar);
	if( (*inchar>=48) && (*inchar<=57) ) { ss>>*obs_state;}		//DLF if ch is 0-9, convert to int
//	if( (*inchar>=48) && (*inchar<=57) ) {*obs_state=atoi(inchar);}		//DLF if ch is 0-9, convert to int
	else if ( (*inchar>=65) && (*inchar<=72) )  {*obs_state=*inchar-55;}//DLF if ch is A-H, convert to 10-17
	else if ( (*inchar>=74) && (*inchar<=78) ) {*obs_state=*inchar-56;}//DLF if ch is J-N, convert to 18-22
	else if ( (*inchar>=80) && (*inchar<=90) ) {*obs_state=*inchar-57;}//DLF if ch is P-Z, convert to 23-33
	
	else if ( (*inchar>=97) && (*inchar<=104) ){*obs_state=*inchar-87;}//DLF if ch is a-h, convert to 10-17
	else if ( (*inchar>=106) && (*inchar<=110) ){*obs_state=*inchar-88;}//DLF if ch is j-n, convert to 18-22
	else if ( (*inchar>=112) && (*inchar<=122) ){*obs_state=*inchar-89;}//DLF if ch is p-z, convert to 23-33
	else return false;
	return true;
}	

void DATA::StateToChr(int *obs_state, char *outchar) {
	stringstream ss;
	
	if (*obs_state<10) {
		ss<<*obs_state;
		strcpy(outchar, ss.str().c_str());
	}
	else {
		if (*obs_state<18) *outchar=*obs_state+55;
		else if (*obs_state<23)  *outchar=*obs_state+56;
		else if (*obs_state<34)  *outchar=*obs_state+57;
		else {
			cout<<"\n\n\t\t****\tError converting state "<< *obs_state <<" to char\n";
			exit(0);
		}
	}
}	

void DATA::MakeNStates() {
	int i, j, k;
	int maxstate=0;

	if (!nchrs) cout << "\n\n***\tCan't get number of states array without nchrs\n";
	else if (!ntaxa) cout << "\n\n***\tCan't get number of states array without ntaxa\n";
	else if (!matrix) cout << "\n\n***\tCan't get number of states array without matrix\n";
	else {
		maxobsstate=0;
		if (!nstates) nstates=new int [*nchrs+1];		
		for (j=1;j<=*nchrs;++j) {
			maxstate=0;
			for (i=1;i<=*ntaxa;++i) {
				for (k=maxstate+1;k<MAXSTATES;++k) if (HasState(i, j, k)) maxstate=k;
			} // for i
			nstates[j]=maxstate+1;		// if the max  state is 1, then there are two states (0 and 1)
			if (chrtype[j]!=4 && nstates[j]>maxobsstate) maxobsstate=nstates[j];	// chrtype is rarely set before this happens in MakeMatrix from file so you gotta call it again somewhere...
		} // for j
	} // else
}

void DATA::MakeStrat() {
	int i, j, k;
	bool fa_is_set;
	
	if (!nchrs) cout << "\n\n***\tCan't get strat ranges array without nchrs\n";
	else if (!ntaxa) cout << "\n\n***\tCan't get strat ranges array without ntaxa\n";
	else if (!matrix) cout << "\n\n***\tCan't get strat ranges array without matrix\n";
	else {
		for (j=1;j<=*nchrs;++j) {					// check for a strat character
			if (chrtype[j]==4) {
				stratchr=j;
				hasStratChr=true;
			}
			if (hasStratChr) break;						// if we've found the strat chr., stop looking (currently only allows one, the first encountered)
		}
		if (hasStratChr) {
			fint = new int [2**ntaxa];
			lint = new int [2**ntaxa];
			for (i=0;i<=*ntaxa;++i) fint[i]=lint[i]=-1;
			for (i=1;i<=*ntaxa;++i)	{				// for each taxon
				fa_is_set=false;
				for (k=0;k<MAXSTATES;++k) {
//					cout << i << "\t" << stratchr << "\t" << matrix[i][stratchr] << endl;
					if (HasState(i, stratchr, k)) {
						if (k>MaxInt) MaxInt=k;
						lint[i]=k;
						if (!fa_is_set) {
							fint[i]=k;
							fa_is_set=true;
						}
					}	
				}
			}
			for (i=*ntaxa+1;i<2**ntaxa;++i) fint[i]=lint[i]=-1;

			for (i=0;i<2**ntaxa;++i) polymorph[i][stratchr]=false;
			hasPolymorph=false;
			for (j=1;j<=*nchrs;++j) 
				if (polymorph[0][j] && chrtype[j]!=4) 
					hasPolymorph=true;	// turns off hasPolymorph if the strat chr is the only polymorphic chr.
		}
	}
}

void DATA::MakeStratSVP() {
	int i, j, k;
	bool fa_is_set, la_is_set;
	
	FAChrPresent=false;
	LAChrPresent=false;
	if (!nchrs) cout << "\n\n***\tCan't get strat ranges array without nchrs\n";
	else if (!ntaxa) cout << "\n\n***\tCan't get strat ranges array without ntaxa\n";
	else if (!matrix) cout << "\n\n***\tCan't get strat ranges array without matrix\n";
	else {
		for (j=1;j<=*nchrs;++j) {					// check for a strat character
			if (chrtype[j]==5) {
				fachr=j;
				FAChrPresent=true;
			}
			else if (chrtype[j]==6) {
				lachr=j;
				LAChrPresent=true;
			}
			if (FAChrPresent && LAChrPresent) break;		// if we've found the strat chr., stop looking (currently only allows one, the first encountered)
		}
		if (!FAChrPresent || !LAChrPresent) {				// if none was found abort
			cout << "\n\n***\tNo FA or LA character found\n";
			return;
		}
		else {
			hasStratChr=true;
			fint = new int [2**ntaxa];
			lint = new int [2**ntaxa];
			for (i=0;i<=*ntaxa;++i) fint[i]=lint[i]=-1;
			for (i=1;i<=*ntaxa;++i)	{				// for each taxon
				fa_is_set=false;
				la_is_set=false;
				for (k=0;k<MAXSTATES;++k) {
					if (!fa_is_set && HasState(i, fachr, k)) {
						if (k>MaxInt) MaxInt=k;
						fint[i]=k;
						fa_is_set=true;
					}	
					if (!la_is_set && HasState(i, lachr, k)) {
						if (k>MaxInt) MaxInt=k;
						lint[i]=k;
						la_is_set=true;			// this would use the youngest LA state if there are more than one state present
					}
					if (fa_is_set && la_is_set) break;
				}
			}
			for (i=*ntaxa+1;i<2**ntaxa;++i) fint[i]=lint[i]=-1;

			for (i=0;i<2**ntaxa;++i) polymorph[i][stratchr]=false;
			hasPolymorph=false;
			for (j=1;j<=*nchrs;++j) 
				if (polymorph[0][j] && chrtype[j]!=4) 
					hasPolymorph=true;	// turns off hasPolymorph if the strat chr is the only polymorphic chr.
		}
	}
}

void DATA::MakeInfChrs() {
	int i, j, k, matches, mismatches, start_taxon;
	
	infchr=new bool[*nchrs+1];
//	infchr[0]=0; 	// this counts the number of inf chrs but only if this is an int array
	start_taxon=1;
	for (j=1;j<=*nchrs;++j) {
		while (start_taxon<*ntaxa && !matrix[start_taxon][j])	++start_taxon;
		i=start_taxon+1;
		matches=mismatches=0;
		while ((matches<2 || mismatches<2) && i<=*ntaxa) {
			for (k=0;k<nstates[j];++k) {
				if ((matrix[start_taxon][j] & matrix[i][j]) & (int) pow(2,k)) ++matches;			// both share the derived state
				else if (!(matrix[start_taxon][j] & (int) pow(2,k)) && (matrix[i][j] & (int) pow(2,k))) ++mismatches;	// state k is observed in taxon j but not start taxon
			}
			++i;
		}
		if (matches<1 || mismatches<2) infchr[j]=false;
		else {
			infchr[j]=true;
//			++infchr[0];	// this counts the number of inf chrs but only if this is an int array
		}
	}
}

void DATA::MakePattern() {
	// MakePattern MUST take place after ChrType and ChrWeights are set
	int i, j, js;
	bool identical;
	
	pattern=new int [*nchrs+1];
	for (j=1;j<=*nchrs;++j) pattern[j]=-1;

	for (j=1;j<*nchrs;++j) {
		if (pattern[j]<0) { 					// if we haven't already considered this chr
			pattern[j]=1;					// consider it a "live" character
			for (js=j+1;js<=*nchrs;++js) {	// for all later characters
				if (pattern[js]<0 && chrtype[j]==chrtype[js] && weight[j]==weight[js]) { // if we haven't already checked the later one, and the character types and weights are identical
					identical=true;				// assume their identical
					for (i=1;i<=*ntaxa;++i) {		// for all taxa
						if ((matrix[i][j] ^ matrix[i][js])!=0) identical=false; // if the states are nonidentical
						if (!identical) break;
					}
					if (identical) {
						pattern[js]=0;		// set the later chr to zero (don't check again)
						++pattern[j];		// increase the number of this chr
					}
				}
			}
		}
	}
	if (pattern[*nchrs]<=0) pattern[*nchrs]=1; 
}

void DATA::PrintMatrix () {
	bool done;
	char *state;

	state=new char [3];
	for (int i=1;i<=*ntaxa;++i) {
		cout.width(3);
		cout << i << " ";
		cout.width(20);
		if (taxon_label) cout<<taxon_label[i]<<"\t";
		else cout << "Taxon " << i << "\t";
		for (int j=1;j<=*nchrs;++j) {
			done=0;
//			if (polymorph[i][j]==1) cout << "{";
			if (isUncertain(i, j)) cout << "{";
			else if (polymorph[i][j]) cout << "(";
			for (int k=0;k<MAXSTATES;++k) {
				if (HasState(i, j, k)) {
					StateToChr(&k, state);
					cout << state;
					done=1;
				}
			}
//			if (polymorph[i][j]==1) cout << "}";
			if (isUncertain(i, j)) cout << "}";
			else if (polymorph[i][j]) cout << ")";
			if (!done) cout <<  "?";
			if (j%10==0) cout << " ";
			}
		cout << "\n";
	}
}

void DATA::PrintFullMatrix () {
	bool done;
	char *state;

	state=new char [3];
	for (int i=1;i<2**ntaxa;++i) {
		cout.width(3);
		cout << i << " ";
		cout.width(20);
		if (i<=*ntaxa) {
			if (taxon_label[i]) cout<<taxon_label[i]<<" ";
			else cout << "Taxon " << i << " ";
		}
		else cout<<"Node "<<i<<" ";
		for (int j=1;j<=*nchrs;++j) {
			done=0;
			if (isUncertain(i, j)) cout << "{";
			else if (polymorph[i][j]) cout << "(";
			for (int k=0;k<*GetMaxObsState();++k) {
				if (HasState(i, j, k)) {
					StateToChr(&k, state);
					cout << state;
					if (matrix[i][j]!=1 && matrix[i][j]%2) cout << " ";
					done=1;
				}
			}
			if (!done) cout <<  "?";
			if (isUncertain(i, j)) cout << "}";
			else if (polymorph[i][j]) cout << ")";
			if (j%10==0) cout << " ";
			}
		cout << "\n";
	}
}

void DATA::PrintSingleChr (int j) {
	bool done;
	char *state;

	state=new char [3];
	for (int i=1;i<2**ntaxa;++i) {
		cout.width(3);
		cout << i << " ";
		cout.width(20);
		if (i<=*ntaxa) {
			if (taxon_label[i]) cout<<taxon_label[i]<<" ";
			else cout << "Taxon " << i << " ";
		}
		else cout<<"Node "<<i<<" ";
		done=0;
		if (isUncertain(i, j)) cout << "{";
		else if (polymorph[i][j]) cout << "(";
		for (int k=0;k<MAXSTATES;++k) {
			if (HasState(i, j, k)) {
				StateToChr(&k, state);
				cout << state;
				done=1;
			}
		}
		if (!done) cout <<  "?";
//		if (polymorph[i][j]==1) cout << "}";
		if (matrix[i][j]>1 && matrix[i][j]%2) cout << "}";
		else if (polymorph[i][j]) cout << ")";
		else if (matrix [i][j]>1 && matrix[i][j]%2!=0) cout << "}";
		cout << "\n";
	}
}

void DATA::PrintStrat () const {
	if (hasStratChr)
		{
		cout << "Taxon\tFA\tLA\n";
		for (int i=1;i<2**ntaxa;++i) {
			if (i<=*ntaxa && taxon_label[i]) cout << taxon_label[i];
			else cout << i;
			cout << "\t\t" << fint[i] << "\t" << lint[i] << "\n";
		}
	}
	else cout << "\n\nNo Stratigraphic Character Present (PrintStrat)\n";
}

void DATA::PrintInclusion () const {
	cout << "Inclusion";
	cout.width(20);
	for (int i=1;i<=*nchrs;++i) {
		cout << included[i];
		if (i%10==0) cout << " ";
	}
	cout << "\n";
}

void DATA::PrintInfChr () const {
	cout << "InfChr";
	cout.width(23);
	for (int i=1;i<=*nchrs;++i) {
		cout << infchr[i];
		if (i%10==0) cout << " ";
	}
	cout << "\n";
}

void DATA::PrintNStates () const {
	cout << "NStates";
	cout.width(23);
	for (int i=1;i<=*nchrs;++i) {
		cout << nstates[i];
		if (i%10==0) cout << " ";
	}
	cout << "\n";
}

void DATA::PrintPolymorph () const {
	int i, j;
	
	cout << "Polymorph:\n";
	for (i=0;i<2**ntaxa;++i) {
		cout.width(3);
		cout << i << " ";
		cout.width(20);
		if (i<=*ntaxa) {
			if (taxon_label[i]) cout<<taxon_label[i]<<" ";
			else cout << "Taxon " << i << " ";
		}
		else cout<<"Node "<<i<<" ";
		for (j=1;j<=*nchrs;++j) {
			cout << polymorph[i][j];
			if (j%10==0) cout << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}

void DATA::PrintNStatesPerTaxon () {
	int i, j, k;
	
	cout << "NStatesPerTaxon:\n";
	for (i=0;i<2**ntaxa;++i) {
		cout.width(3);
		cout << i << " ";
		cout.width(20);
		if (i<=*ntaxa) {
			if (taxon_label[i]) cout<<taxon_label[i]<<" ";
			else cout << "Taxon " << i << " ";
		}
		else cout<<"Node "<<i<<" ";
		for (j=1;j<=*nchrs;++j) {
			k=(float) getNumStates(i,j);
			cout << k;
			if (j%10==0) cout << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}

bool DATA::HasState(int taxon, int chr, int req_state) {
	unsigned long long x;
	
	x=(unsigned long long) (pow(2, req_state));
	if (matrix[taxon][chr] & x) return true;
	else return false;
}

bool DATA::ShareAnyState(int taxon1, int taxon2, int chr) {
	if (matrix[taxon1][chr] & matrix[taxon2][chr]) return true; 
	else return false;
}

bool DATA::IsMissing(int taxon1, int chr) {
	if (matrix[taxon1][chr]) return false;
	else return true;
}

bool DATA::isUncertain(int taxon, int chr) {
	int i;
	bool flag=false;
	
	if (matrix[taxon][chr] && !polymorph[taxon][chr]) {	// if the chr is not missing nor polymorphic
		for (i=0;i<nstates[chr];++i) {
			if (matrix[taxon][chr] & (unsigned long long) (pow(2, i))) {
				if (!flag) flag=true;
				else return true;
			}
		}
	}
	return false;
}

int DATA::getNumStates(int taxon, int chr) {
	int i, num_states=0;
	for (i=0;i<nstates[chr];++i) if (matrix[taxon][chr] & (unsigned long long) (pow(2, i))) ++num_states;
	return num_states;
}

int DATA::MinState(int taxon, int chr) {
	for (int i=0;i<nstates[chr];++i) {
		if (matrix[taxon][chr] & (unsigned long long) (pow(2, i))) return i;
	}
	return 0;
}

int DATA::MaxState(int taxon, int chr) {
	for (int i=nstates[chr]-1;i>=0;--i) { 
		if (matrix[taxon][chr] & (unsigned long long) (pow(2, i))) return i;
	}
	return 0;	
}

int DATA::MinSharedState(int taxon1, int taxon2, int chr) {
	int i;
	unsigned long long x;
	bool state_found=false;
	
	x = matrix[taxon1][chr] & matrix[taxon2][chr];
	for (i=0;i<nstates[chr];++i) {
		if (x & (unsigned long long) (pow(2, i))) state_found=true;
		if (state_found) return i;
	}
	return 0;
}

void DATA::BootstrapMatrixByTaxa()
{
	int i, j;
	float x;
	int selected_sp;
	
	if (!boot_matrix) boot_matrix=ullmatrix(*ntaxa+1, *nchrs+1);
	for (i=1;i<=*ntaxa;++i) {
		x=0;
		x=((float)rand()/(float)RAND_MAX)*(*ntaxa);
		selected_sp=(int)x+1;
		for (j=1;j<=*nchrs;++j) boot_matrix[i][j]=matrix[selected_sp][j];
	}
}

void DATA::BootstrapMatrixByChrs() {
	int i, j;
	float x;
	int selected_chr;
	
	if (!boot_matrix) boot_matrix=ullmatrix(*ntaxa+1, *nchrs+1);
	for (j=1;j<=*nchrs;++j) {
		x=0;
		x=((float)rand()/(float)RAND_MAX)*(*nchrs);
		selected_chr=(int)x+1;
		for (i=1;i<=*nchrs;++i) boot_matrix[i][j]=matrix[i][selected_chr];
	}
}

int ** DATA::outputStratAsMatrix() {
	int **strat=imatrix(*ntaxa, 2);
	for (int i=1;i<=*ntaxa;i++) {
		strat[i-1][0]=*GetFInt(i);
		strat[i-1][1]=*GetLInt(i);
	}
	return strat;
}

