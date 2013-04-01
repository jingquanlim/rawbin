#ifndef __EXTEND_GUARD__
#define __EXTEND_GUARD__

#include "file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct Junction {
	unsigned p;
	unsigned q;
	unsigned r;
	char *Chrom;
	int Mismatches;
	char signal[5];
	int isCanonical();
	int Label;
	int ID;
	char Sign;
	float score;
};
extern unsigned char* Original_Text;
int Canonical_Score(char* signal);
Junction* extend(char* R, unsigned x, unsigned y, unsigned p, unsigned q, char sign);
void loadPac(char* filename);
void Get_Bases_ASCII (unsigned Location,int StringLength,char* Org_String);
void Get_Bases(unsigned Location,int StringLength,char* Org_String);
inline void Convert_Reverse(char* Read_bin, char* RC_bin,int StringLength);
float signalScore(char* signal);
float getScore(unsigned p, unsigned q, unsigned x, int size, int i, int misL, int misR, char sign);
#endif
