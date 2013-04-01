//BATMAN 1.1 - added Substring search..
//BATMAN 1.10 - enhanced default output processing...
//		handle N's.
//		print blanks...
//		--maxhits bug fixes
//{-----------------------------  INCLUDE FILES  -------------------------------------------------/
#include <stdio.h>
#include <string>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include "bfix.h"
#include <getopt.h>
#include "zlib.h"
#include "const.h"
#include "assert.h"
#include "Hash.h"
#include "rqindex.h"
#include <fstream>
#include <iostream>
#include "Print.h"
#include "init.h"
#include <vector> 
#include <math.h>
extern "C" 
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}
//}-----------------------------  INCLUDE FILES  -------------------------------------------------/

#include "batlib.h"
#include "Cmdline.h"
#include "Indexes.h"
#include "file.h"
#include "extend.h"
//#include "batlibX.cpp"
//#include "rqindex.cpp"
//{-----------------------------  DEFINES  -------------------------------------------------/

//
//{---------------------------- GLOBAL VARIABLES -------------------------------------------------
extern const int UNIQUE_SIGNAL=0;
extern const int UNIQUE_NOSIGNAL=1;
extern const int NON_UNIQUE_SIGNAL=2;
extern const int NON_UNIQUE_OTHERS=3;
extern const int MINX=18;//Minimum extension..
extern const int MIS_DENSITY=1;

const int DUMMY_JUNC=1000;
const int TRANCRIPT_END=1001;
extern const bool DEBUG=false;

struct Transcript
{
	unsigned x;
	unsigned y;
};

//Transcript Trans_Array[20];
const int MAX_JUNCS_ALLOWED=1000;
const int MAX_JUNCS_IN_TRANSCRIPT=20;//Maximum number of junctions that are allowed in a trascript..
const int MAX_INSPECTED_PAIRS=20000;//INT_MAX;
Junction Trans_Array[MAX_JUNCS_IN_TRANSCRIPT];
Junction Compiled_Junctions[MAX_JUNCS_ALLOWED];
int Compiled_Junctions_Ptr;

int INDEX_RESOLUTION=30000;
int EXONGAP;
int READLEN;
unsigned CONVERSION_FACTOR;
char* LOG_SUCCESS_FILE=NULL;
FILE *Log_SFile;
int MAX_MISMATCHES=1;
BWT *fwfmi,*revfmi;
unsigned Conversion_Factor;
map <unsigned, Ann_Info> Annotations;
unsigned Location_Array[80];
unsigned char* Original_Text;
char Char_To_CodeC[256];
char Char_To_Code[256];
Offset_Record Genome_Offsets[80];
unsigned Offsets[80];
int ORG_STRINGLENGTH;

unsigned RQ_Hits;
SA* SA_Index;
char COMPRESS;
int *SA_Blocks;
int Max_Junc_Found;
int Least_Mis_In_Junc;
MEMX Generic_Hits;

float Donor_Prob[16][64][2];
float Acc_Prob[16][64][2];
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------



//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*
int Find_All_Single_Junc(char* Read,char* Converted_Read,MEMX & Pre,MEMX & Suf,int STRINGLENGTH,LEN & L,PAIR* Pairs,Junction *Final_Juncs,int & Err,int & Label);
bool Align(char* Source,int StringLength,unsigned Dest_Loc,SARange & R,int Actual_Length,int Read_Skip);
int Seek_Junc(char* S,SARange R,int Read_Skip,int Junc_Count,int Mis_In_Junc_Count,unsigned Last_Exon,int StringLength,int Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int & Inspected_Pairs,int & Err,bool Sign);
int Seek_Single_Strand(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,int Sign);
int Seek_All_Junc(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf);
int Find_Single_Junc(char* Read,char* Read_Head,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err,int & Label,char Sign,unsigned Anchor,SARange & R,int Skip);
void Enum_Single_Junctions(char* Org_Read,char* Converted_Read,int Read_Skip,int StringLength, unsigned Anchor,int & Inspected_Pairs,SARange & R,int Junc_Count,int Mis_In_Junc_Count,int & Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int Skip,int & Err,int Sign);
void Load_All_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,RANGEINDEX & Range_Index);
void Init(BWT *revfmi,unsigned & SOURCELENGTH,PAIR* & Pairs,gzFile & Input_File,gzFile & Mate_File,FILETYPE & File_Info,Parameters & CL,FILE* & OUT,Index_Info & Genome_Files);
bool  Progress_Bar(Parameters & CL,unsigned & Number_of_Tags,unsigned & Progress,unsigned & Tag_Count,FILETYPE & File_Info);
int Find_Pairings(int & Pairs_Index,SARange* MF_Pre_Hits,SARange* MF_Suf_Hits,PAIR* &  Pairs, char* Read,char* Read_Head,int STRINGLENGTH,int & Final_Juncs_Ptr,Junction *Final_Juncs,int & Min_Mismatch,int Mis_In_Head,int Mis_In_Tail,int Skip,int & Label,char Sign);
void Open_Outputs(ofstream & SAM,string filename);
inline char* Nullify_String(char* S);
void Print_Hits(READ & Head,Junction *Final_Juncs,FILE* OUT,ofstream & SAM,int Tag_Count, int firstSignal);
int Classify_Hits(Junction * Final_Juncs, int & firstSignal);
inline void Convert_Reverse(char* Read,char * RC_Read,char* RC_Bin,int StringLength);
void Print_SAM_Header(std::map <unsigned, Ann_Info> Annotations,int argc,char* argv[],char* Input_File);
void Init_Prob();
void fillProbArray(float prob[16][64][2], string filename);
int getBest(char* Current_Tag,int StringLength,Junction * Compiled_Junctions, int * approved, bool print);
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

int main(int argc, char* argv[])
{
	unsigned Total_Hits=0,Tags_Processed=0,Tag_Count=0;
	time_t Start_Time,End_Time;
	unsigned Number_of_Tags=1000;
	unsigned Progress=0;
	FILE* OUT;

//------------------- INIT -------------------------------------------------------------
	Index_Info Genome_Files;
	Parameters CL;	
	MMPool *mmPool;
	RANGEINDEX Range_Index;
	unsigned SOURCELENGTH;
	ofstream SAM[4];
	string SAMFiles[] = {"unique_signal.sam","unique_nosignal.sam",/*"mishits.fq",*/"nonunique_onesignal.sam","others.sam"};
	PAIR *Pairs;
	//Junction Final_Juncs[2*MAX_JUNCS_TO_STORE+2];
	gzFile Input_File,Mate_File;
	FILETYPE File_Info;

	Parse_Command_line(argc,argv,Genome_Files,CL);
	Load_All_Indexes(Genome_Files,fwfmi,revfmi,mmPool,Range_Index);
	Init(revfmi,SOURCELENGTH,Pairs,Input_File,Mate_File,File_Info,CL,OUT,Genome_Files);
	Open_Genome_Files(Genome_Files.LOCATIONFILE,Genome_Offsets,Offsets);
	for(int i=0;i<4;i++)
	{
		Open_Outputs(SAM[i],SAMFiles[i]);
	}
	ofstream rejectedSAM;
	string rejectedSAMFile = "rejected.sam";
	Open_Outputs(rejectedSAM, rejectedSAMFile);
	
	Print_SAM_Header(Annotations,argc,argv,CL.PATTERNFILE);
//------------------- INIT -------------------------------------------------------------

	time(&Start_Time);
	READ Head,Tail;
	int LOOKUPSIZE=3;
	MEMLOOK MLook;MLook.Lookupsize=3;
	Build_Tables(fwfmi,revfmi,MLook);
	LEN L;L.IGNOREHEAD=0;
	Split_Read(RQFACTOR,L);//we are scanning 18-mers...
	Conversion_Factor=revfmi->textLength-RQFACTOR;
//--------------------- Setup Data Structure for Batman Prefix ----------------------------------------
	MEMX MF_Pre,MF_Suf;//MemX is the data structure for doing Batman alignment. MF_Pre is for the prefix, MC for suffix..
	Init_Batman(MF_Pre,L,MLook,MAX_MISMATCHES);
	Init_Batman(MF_Suf,L,MLook,MAX_MISMATCHES);
	Init_Batman(Generic_Hits,L,MLook,MAX_MISMATCHES);
//--------------------- Setup Data Structure for Batman End----------------------------------------
	

//--------------------- Load probability information --------------------------
	Init_Prob();

	time_t Maptime;
	{
		fprintf(stderr,"======================]\r[");//progress bar....
		int Actual_Tag=0;
		int selectedJunctions[MAX_JUNCS_ALLOWED];
		while (Read_Tag(Head,Tail,Input_File,Mate_File,File_Info))
		{
			if(!Progress_Bar(CL,Number_of_Tags,Progress,Tag_Count,File_Info)) break;
			if(CL.MAX_TAGS_TO_PROCESS && CL.MAX_TAGS_TO_PROCESS<Actual_Tag) break;

			int Label=0;
			Actual_Tag++;

			Max_Junc_Found=INT_MAX;
			Least_Mis_In_Junc=INT_MAX;
			Compiled_Junctions_Ptr=0;
			int Err=Seek_All_Junc(Head.Tag,File_Info.STRINGLENGTH,MF_Pre,MF_Suf);
			Compiled_Junctions[Compiled_Junctions_Ptr].p=UINT_MAX;

			if(Compiled_Junctions_Ptr && !(DEBUG && Err))
			{
				int firstSignal = -2;
				int tempType = Classify_Hits(Compiled_Junctions,firstSignal);
				int approvedPtr;
				approvedPtr = getBest(Head.Tag,File_Info.STRINGLENGTH,Compiled_Junctions, selectedJunctions, true);
				
				if(approvedPtr > 1) {
					for(int i=0; i<approvedPtr; i++) {
						Print_Hits(Head,Compiled_Junctions,OUT,rejectedSAM,Tag_Count,selectedJunctions[i],tempType);//tempType);
					}
				}	
				else {
					for(int i=0; i<approvedPtr; i++) {
						Print_Hits(Head,Compiled_Junctions,OUT,SAM[tempType],Tag_Count,selectedJunctions[i],tempType);//tempType);
					}
				}

			}
			else
			{
				//if(Err) SAM[1] << ">ERR:";
				//SAM[1] << Head.Description << Head.Tag_Copy;

			}

			/*if(Find_All_Single_Junc(Head.Tag_Copy,Head.Tag,MF_Pre,MF_Suf,File_Info.STRINGLENGTH,L,Pairs,Final_Juncs,Err,Label))
			{
				if(Label==1)//One junction..
				{
					int firstSignal = -2;
					int tempType = Classify_Hits(Head,Final_Juncs,firstSignal);
					if(tempType != 2)
						firstSignal = -2;
					Print_Hits(Head,Final_Juncs,OUT,SAM[tempType],Tag_Count,firstSignal,tempType);
				}
				else//Multiple Junctions..
				{
					
				}
			}
			else//Read is not mapped
			{
				
				if(!Err) ;//printf(">%d-%s%s",Tag_Count,Head.Description,Head.Tag_Copy);
			}*/
		}
		UnLoad_Indexes(fwfmi,revfmi,mmPool,Range_Index);
		fprintf(stderr,"\r[++++++++100%%+++++++++]\n");//progress bar....
		time(&End_Time);
		Maptime=difftime(End_Time,Start_Time);

		char Junction_File[]="junctions.bed";
		Print_Junctions(Junction_File);

	}
	/*printf("Average %d\n",Sum);*/
	fprintf(stderr,"%u / %u Tags/Hits",Tag_Count,Total_Hits);
	time(&End_Time);fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",(float)difftime(End_Time,Start_Time));
	fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",float(Maptime));
	return 0;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Convert Reverse
 *  Description:  Convert Read in binary form to its binary and nucleotide rev. comp.
 * =====================================================================================
 */

inline void Convert_Reverse(char* Read,char * RC_Read,char* RC_Bin,int StringLength)
{
	for (unsigned i=0;i<=StringLength-1;i++)
	{
		RC_Bin[StringLength-1-i]=3-Read[i];
		RC_Read[StringLength-1-i]="ACGT"[3-Read[i]];
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Find_All_Single_Junc
 *  Description:  Find all Single junctions in the read Read and its reverse compliment.. 
 *  		  Returns least mismatch junctions in Final_Juncs, terminated by sentinel Final_Junc[Last].p=INT_MAX 
 *  		  Err is set if there were overflows in hits, 
 *  		  returns zero if no junctions found.. 
 *  		  Label will increase by one for each putative junction, disregarding those found in extensions
 * =====================================================================================
 */
/*int Find_All_Single_Junc(char* Read,char* Converted_Read,MEMX & Pre,MEMX & Suf,int STRINGLENGTH,LEN & L,PAIR* Pairs,Junction *Final_Juncs,int & Err,int & Label)
{
	//int Result=Find_Single_Junc(Head.Tag_Copy,Head.Tag,MF_Pre,MF_Suf,File_Info.STRINGLENGTH,L,Pairs,Final_Juncs,Err,Label);
	char RC_Read[MAXDES],RC_Bin[MAXDES];
	int Result1=Find_Single_Junc(Read,Converted_Read,Pre,Suf,STRINGLENGTH,L,Pairs,Final_Juncs,Err,Label,'+');
	Convert_Reverse(Converted_Read,RC_Read,RC_Bin,STRINGLENGTH);Final_Juncs+=Result1;
	int Result2=Find_Single_Junc(RC_Read,RC_Bin,Pre,Suf,STRINGLENGTH,L,Pairs,Final_Juncs,Err,Label,'-');
	return Result1+Result2;
}*/


int Classify_Hits(Junction * Final_Juncs, int & firstSignal){
	int signalCount = 0, count = 0;
	for(int i=0; Final_Juncs[i].p != UINT_MAX; i++){
		count++;
		if(Final_Juncs[i].q!=Final_Juncs[i].p-1)//not an exact match?
		{
			if(Final_Juncs[i].isCanonical()) {
				signalCount++;
				if(firstSignal < 0)
					firstSignal = i;
			}
		}
	}
	//printf("Classify count: %d\n", count);
	if(count == 1 && signalCount == 1)
		return UNIQUE_SIGNAL;
	else if(count == 1)
		return UNIQUE_NOSIGNAL;
	else if(count > 1 && signalCount == 1)
		return NON_UNIQUE_SIGNAL;
	else
		return NON_UNIQUE_OTHERS;
}

int getBest(char* Current_Tag,int StringLength,Junction * Final_Juncs, int * approved, bool print) {
	float max = -100000;
	int ptr = 0;
	int count = 0;
	char RC_Read[MAXDES],RC_Bin[MAXDES];
	char Cat[MAXDES];
	Convert_Reverse(Current_Tag,RC_Read,RC_Bin,StringLength);

//Calculate true mismatches..
	for(int i =0; Final_Juncs[i].p != UINT_MAX; i++) 
	{
		if(Final_Juncs[i].q)
		{
			Get_Bases_ASCII(Final_Juncs[i].p-Final_Juncs[i].r, Final_Juncs[i].r, Cat);
			Get_Bases_ASCII(Final_Juncs[i].q+1, StringLength-Final_Juncs[i].r, Cat+Final_Juncs[i].r);
		}
		else
		{
			Get_Bases_ASCII(Final_Juncs[i].p, StringLength,Cat);
		}

		char* Read;
		Final_Juncs[i].Mismatches=0;
		if(Final_Juncs[i].Sign)
			Read=Current_Tag;
		else
			Read=RC_Bin;

		for(int j=0;j<StringLength;j++)
		{
			if(Cat[j]!="ACGT"[Read[j]]) 
				Final_Juncs[i].Mismatches++;
			Cat[j]=0;
		}
	}

	for(int i =0; Final_Juncs[i].p != UINT_MAX; i++) {
		count ++;
		if(!Final_Juncs[i].q)
			Final_Juncs[i].score = 5;
		float tempScore = Final_Juncs[i].score - Final_Juncs[i].Mismatches;
		if(fabs(tempScore - max) < 0.0000001) {
			approved[ptr] = i;
			++ptr;
		}
		else if(tempScore > max) {
			max = tempScore;
			ptr = 0;
			approved[ptr] = i;
			++ptr;
		}	
			//if(print)
				//printf("%f\t%u\t%u\t%u\n",Final_Juncs[i].score,Final_Juncs[i].p, Final_Juncs[i].q, Final_Juncs[i].r);
	}
	//printf("Total reported: %d\n",count);
	return ptr;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Find_Single_Junc
 *  Description:  Find Single junctions in the read Read 
 *  		  if Anchor is valid, it gives left anchor of an 18-mer..
 *  		  Returns least mismatch junctions in Final_Juncs, terminated by sentinel Final_Junc[Last].p=UINT_MAX 
 *  		  Err is set if there were overflows in hits, 
 *  		  returns zero if no junctions found.. 
 *  		  Label will increase by one for each putative junction, disregarding those found in extensions
 * =====================================================================================
 */
int Find_Single_Junc(char* Read,char* Read_Head,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err,int & Label,char Sign,unsigned Anchor,SARange & R,int Skip)
{
	assert(STRINGLENGTH<=READLEN && Read && Converted_Read && Skip>=0 && STRINGLENGTH >0 && (Sign=='+'||Sign=='-'));
	int Pairs_Index=0;
	int Final_Juncs_Ptr=0;
	int Ret_Value=0;//if Any err in pairing, gets set..
	int Min_Mismatch=INT_MAX;
	//Get Hits for head..
	MF_Pre.Hits=0;MF_Pre.Hit_Array_Ptr=0;MF_Pre.Current_Tag=Converted_Read;MF_Pre.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
	int Last_Mis_Suf= -1;
	int Last_Mis_Pre=0;
	int Mis_In_Head=R.Mismatches;
	if(Anchor==UINT_MAX)//If not extending from a putative exon..
	{
		if(Skip==0)// || (R.End-R.Start < SAGAP_CUTOFF_T))//first l-mer..
		{
			MF_Pre.Hit_Array[0]=R;MF_Pre.Hit_Array_Ptr=1;//Read_Head=0;
		}
		else
		{
			MF_Pre.Current_Tag=Converted_Read;
			Last_Mis_Pre=Scan(MF_Pre,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
		}
	}
	else
	{
		MF_Pre.Hit_Array[0].Start=MF_Pre.Hit_Array[0].End=Anchor;MF_Pre.Hit_Array_Ptr=1;
		Mis_In_Head=R.Mismatches;//Read_Head=0;
	}
	int Top_Hit_End=MF_Pre.Hit_Array_Ptr;//Save the end point top hits..
	//------------------------------------ starts Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	if(Last_Mis_Pre>=0)//Did we get a hit?
	{
		//Get Hits for suffix..
		MF_Suf.Hits=0;MF_Suf.Hit_Array_Ptr=0;MF_Suf.Current_Tag=Converted_Read+(STRINGLENGTH-RQFACTOR);MF_Suf.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
		Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
		int Top_Hit_End_Suf=MF_Suf.Hit_Array_Ptr;//Save the end point top hits..
	}
	if(Last_Mis_Pre >=0 && Last_Mis_Suf>=0)
	{
		Err+=Find_Pairings(Pairs_Index,MF_Pre.Hit_Array,MF_Suf.Hit_Array,Pairs,Read,Read_Head,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,Mis_In_Head,Last_Mis_Suf,Skip,Label,Sign);
	}
	else return 0;//no anchors
	/*if (Final_Juncs_Ptr || Err)
	{
		Final_Juncs[Final_Juncs_Ptr].p=UINT_MAX;
		return Final_Juncs_Ptr;
	}*/
	//------------------------------------ End Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	//------------------------------------ Start Higest/Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	SARange* Sub_Opt_Pre_Start=MF_Pre.Hit_Array+MF_Pre.Hit_Array_Ptr;
	SARange* Sub_Opt_Suf_Start=MF_Suf.Hit_Array+MF_Suf.Hit_Array_Ptr;
	int M1=Mis_In_Head,M2=Last_Mis_Suf;

	Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis_Suf+1,UINT_MAX);
	if(Anchor==UINT_MAX)
	{
		Last_Mis_Pre=Scan(MF_Pre,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis_Pre+1,UINT_MAX);
	}
	else
	{
		Last_Mis_Pre= -1;//fake no subopt hit..
	}
	if(Last_Mis_Pre >=0)
	{
		Err+=Find_Pairings(Pairs_Index,Sub_Opt_Pre_Start,MF_Suf.Hit_Array,Pairs,Read,Read_Head,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,Mis_In_Head,M2,Skip,Label,Sign);
	}
	if(Last_Mis_Suf >=0)
	{
		Err+=Find_Pairings(Pairs_Index,MF_Pre.Hit_Array,Sub_Opt_Suf_Start,Pairs,Read,Read_Head,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,Mis_In_Head,Last_Mis_Suf,Skip,Label,Sign);
	}
	//------------------------------------ End Higest/Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	//------------------------------------ End Higest Mismatch Pairing of Suf/Pref -------------------------------------------------
	if(Last_Mis_Pre >=0 && Last_Mis_Suf>=0)
	{
		Err+=Find_Pairings(Pairs_Index,Sub_Opt_Pre_Start,Sub_Opt_Suf_Start,Pairs,Read,Read_Head,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,Mis_In_Head,Last_Mis_Suf,Skip,Label,Sign);
	}
	//------------------------------------ End Higest Mismatch Pairing of Suf/Pref -------------------------------------------------

	Final_Juncs[Final_Juncs_Ptr].p=UINT_MAX;
	//if (Err) return 0;
	return Final_Juncs_Ptr;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Find_Pairings
 *  Description:  Find Junctions by pairing SA-Ranges in MF_Pre and MF_Suf.. 
 *  		  Returns least mismatch junctions, with mismatches not exceeding Min_Mismatch in Final_Juncs. Last junc at Final_Junc[Final_Juncs_Ptr-1] 
 *  		  Err is set if there were overflows in hits, and returns zero..
 * =====================================================================================
 */
struct Split_Map
{
	int p;
	int q;
	int r;
};

int Find_Pairings(int & Pairs_Index,SARange* MF_Pre_Hits,SARange* MF_Suf_Hits,PAIR* &  Pairs, char* Read,char* Read_Head,int STRINGLENGTH,int & Final_Juncs_Ptr,Junction *Final_Juncs,int & Min_Mismatch,int Mis_In_Head,int Mis_In_Tail,int Skip,int & Label,char Sign)
{
	assert(STRINGLENGTH<=READLEN && Skip>=0 && STRINGLENGTH >0 && (Sign=='+'||Sign=='-'));

	int Ret_Value=0;//if Any err in pairing, gets set..
	Pairs_Index=0;
	int Err=Pair_Reads(MF_Pre_Hits,MF_Suf_Hits,Pairs,Pairs_Index);
	if(!Pairs_Index) return Err;	

	while(--Pairs_Index>=0) 
	{
		Ann_Info A,A1;
		unsigned L1=Pairs[Pairs_Index].Head;
		unsigned L2=Pairs[Pairs_Index].Tail;					
		assert (L1>=0 && L2>=0);
		if(L2-L1<STRINGLENGTH-RQFACTOR) continue; //Too close.. can be an indel though..

		if(Skip && Read_Head)//Check if the beginning of the strings matches..
		{
			char Str[200];
			Get_Bases_ASCII(L1-Skip,Skip,Str);
			bool Head_Match=true;int Mis_Count=0;
			for(int i=0;Read_Head[i];i++)
			{
				assert(Read_Head[i]>='A' && Read_Head[i]<='T' && Str[i]>='A' && Str[i]<='T');
				if(Read_Head[i]!=Str[i])
				{
					Mis_Count++;
					if(Mis_Count>Mis_In_Head)
					{
					       	Head_Match=false;break;
					}
				}
			}
			if(!Head_Match) continue;
		}

		Junction* junctions = extend(Read,MINX-1,STRINGLENGTH-MINX,L1,L2,Sign);
		if(junctions[0].p==UINT_MAX) continue;

		Location_To_Genome(Pairs[Pairs_Index].Head,A);
		Location_To_Genome(Pairs[Pairs_Index].Tail,A1);
		if (Pairs[Pairs_Index].Head+STRINGLENGTH > A.Size||Pairs[Pairs_Index].Tail+STRINGLENGTH > A.Size)//check for a Boundary Hit..
		{
			continue;
		}
		if (A1.ID == A.ID)//in the same chromosome?
		{
			if(DEBUG || junctions[0].Mismatches<=Min_Mismatch)
			{
				if(!DEBUG && Min_Mismatch>junctions[0].Mismatches) 
				{
					Min_Mismatch=junctions[0].Mismatches;Final_Juncs_Ptr=0;//Store min mismatch junctions...
				}
				int copyCount = 0;
				for(int i=0;junctions[i].p!=UINT_MAX;i++)
				{
					copyCount ++;
					Split_Map S;
					if(Final_Juncs_Ptr>=MAX_JUNCS_TO_STORE-1) {Err=1;break;}
					//Location_To_Genome(junctions[i].p,A);
					//Location_To_Genome(junctions[i].q,A);
					junctions[i].r+=Skip;
					S.p=junctions[i].p;
					S.q=junctions[i].q; 	
					S.r=junctions[i].r;

					junctions[i].Chrom=A.Name;
					junctions[i].Label=Label;
					junctions[i].Sign=Sign;
					junctions[i].ID=A.ID;

					FILE* F=Genome_Offsets[A.ID].Out_File;
					fwrite(&S,sizeof(Split_Map),1,F);

					//if(S.q!=S.p-1)//not an exact match?
					{
						Final_Juncs[Final_Juncs_Ptr++]=junctions[i];
					}
				}
				//printf("copy count: %d\n", copyCount);
			}
		}
		//printf("Final_Junc_Ptr: %d\n",Final_Juncs_Ptr);
		delete [] junctions;
	}
	Label++;
	return Err;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Align
 *  Description:  Matches S[Source,StringLength] to G[Dest,...] allowing k mis (if Dest!=UINT_MAX) else find S[Source,StringLength] in G allow
 *  		  1 mis as SARanges.. 
 *		  SARanges reported in Generic_Hits
 *  		  Returns true if found.. 
 * =====================================================================================
 */

bool Align(char* Source,int StringLength,unsigned Dest_Loc,SARange & R,int Actual_Length,int Read_Skip)
{
	if(Dest_Loc!=UINT_MAX)//extend in genome..
	{
		char Dest[StringLength];
		Get_Bases(Dest_Loc,StringLength,Dest);
		int Mis_Count=0;
		for(int i=0;i<StringLength && Mis_Count<=1/*Max mis allowed*/;i++)
		{
			assert(Dest[i]<=4 && Dest[i]>=0 && Source[i]<=4 && Source[i]>=0);
			if(Dest[i]!=Source[i])
			{
				Mis_Count++;
			}
		}
		if (Mis_Count>MIS_DENSITY)
			return false;
		else
		{
			R.Mismatches+=Mis_Count;
			Generic_Hits.Hit_Array[0]=R;
			Generic_Hits.Hit_Array_Ptr=1;
			return true;
		}
	}	
	else
	{
		Generic_Hits.Hit_Array_Ptr=0;R.Level=1;
		int T_Len=Generic_Hits.L.STRINGLENGTH;Generic_Hits.L.STRINGLENGTH=Actual_Length;
		int Lo=0;
		Extend_Forwards(Source,R,/*MAX_MISMATCHES*/MIS_DENSITY,1,StringLength,INT_MAX,revfmi,Generic_Hits,true,Read_Skip,Lo);
		if(Lo)
		{
			SARange SA;
			SA=Generic_Hits.Hit_Array[Lo];
			Generic_Hits.Hit_Array[Lo]=Generic_Hits.Hit_Array[0];
			Generic_Hits.Hit_Array[0]=SA;
			assert(SA.Mismatches==R.Mismatches);
		}
		Generic_Hits.L.STRINGLENGTH=T_Len;
		if(Generic_Hits.Hit_Array_Ptr) return true;
		else return false;
	}
}

inline void B2C(char *S,char *D,int StringLength)
{
	*(D+StringLength)=0;
	for(int i=0;i<StringLength;i++,S++,D++)
		*D="ACGT"[*S];
}

int Seek_Junc(char* S,SARange R,int Read_Skip,int Junc_Count,int Mis_In_Junc_Count,unsigned Last_Exon,int StringLength,int Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int & Inspected_Pairs,int & Err,bool Sign)
{
	assert(Trans_Array_Ptr<MAX_JUNCS_IN_TRANSCRIPT);assert(Read_Skip>=0 && StringLength >=0 && StringLength<=READLEN && Inspected_Pairs>=0);
	if(Read_Skip>StringLength-MINX)//End of read..
	{
		if(Junc_Count) 
		{
			Max_Junc_Found=Junc_Count;
		}
		if(Align(S+Read_Skip,StringLength-Read_Skip,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip,R,StringLength,Read_Skip))//End of the read can be matched..
		{
			if(!Junc_Count)//Match..
			{
				Junction J;J.Sign=(Sign)? 1:0;J.Mismatches=R.Mismatches;
				if (Last_Exon!=UINT_MAX) 
				{
					J.p=Last_Exon;
					J.r=J.q=0;
					if(Compiled_Junctions_Ptr<MAX_JUNCS_ALLOWED)
					{
						Compiled_Junctions[Compiled_Junctions_Ptr++]=J;
					}
					else 
					{
						Err=1;
					}
				}
				else
				{
					for(unsigned i=R.Start;i<=R.End;i++)
					{
						J.p=revfmi->textLength-READLEN-BWTSaValue(revfmi,i);;
						J.r=J.q=0;
						if(Compiled_Junctions_Ptr<MAX_JUNCS_ALLOWED)
						{
							Compiled_Junctions[Compiled_Junctions_Ptr++]=J;
						}
						else 
						{
							Err=1;break;
						}
					}
				}
				return TRANCRIPT_END;
			}

			if(Least_Mis_In_Junc> Mis_In_Junc_Count) Least_Mis_In_Junc=Mis_In_Junc_Count; 
			for(int i=0;i<Junc_Count;i++)
			{
				if(Compiled_Junctions_Ptr< MAX_JUNCS_ALLOWED)
				{
					Trans_Array[i].Sign=(Sign)? 1:0;Trans_Array[i].Mismatches=R.Mismatches;Compiled_Junctions[Compiled_Junctions_Ptr++]=Trans_Array[i];
				}
				else
				{
					Err=1;
				}
			}
			return TRANCRIPT_END;
		}
		else//End of the read may contain extra junc that cannot be processed..
		{
			if(Junc_Count && Least_Mis_In_Junc> Mis_In_Junc_Count) Least_Mis_In_Junc=Mis_In_Junc_Count; 
			for(int i=0;i<Junc_Count;i++)
			{
				if(Compiled_Junctions_Ptr< MAX_JUNCS_ALLOWED)
				{
					Trans_Array[i].Sign=(Sign)? 1:0;Trans_Array[i].Mismatches=R.Mismatches;Compiled_Junctions[Compiled_Junctions_Ptr++]=Trans_Array[i];
				}
				else
				{
					Err=1;
				}
			}
			return DUMMY_JUNC;
		}
	}

	if(Align(S+Read_Skip,MINX,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip,R,Read_Skip+MINX,Read_Skip))//Can exact extension be done..
	{
		std::vector <SARange> Hits(Generic_Hits.Hit_Array_Ptr);
		for (int i=0;i<Generic_Hits.Hit_Array_Ptr;i++)//Save hits in a vector..
		{
			assert(Generic_Hits.Hit_Array[i].Start);
			Hits[i]=Generic_Hits.Hit_Array[i];
		}
		unsigned Last_ExonT;int Hit_Count=Generic_Hits.Hit_Array_Ptr-1;
		for (int i=0;i<=Hit_Count;i++)//Save hits in a vector..
		{
			SARange SA=Hits[i];
			Last_ExonT=Last_Exon;
			if(SA.Start==SA.End && Last_Exon==UINT_MAX) 
			{
				Last_ExonT=SA.Start-Read_Skip;
			}
			Seek_Junc(S,SA,Read_Skip+MINX,Junc_Count,Mis_In_Junc_Count,Last_ExonT,StringLength,Trans_Array_Ptr,Pre,Suf,Inspected_Pairs,Err,Sign);
			if (Inspected_Pairs >= MAX_INSPECTED_PAIRS) {Err++;return 0;}
		}
	}

//Junction check routines
	if(Read_Skip>=MINX)
	{
		//Case of a single junction between S[a,a+2l-1]
		if (Junc_Count<Max_Junc_Found) //This will add one junction..
		{
			//the junction is between S[Read_Skip,Read_Skip+18]
			//so anchor S[Read_Skip-18,Read_Skip]
			//Penguin S[Read_Skip-18,Read_Skip] ---- S[Read_Skip+18,Read_Skip+36] and extend..
			Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign);
		}
	}

	return Trans_Array_Ptr;

}

void Enum_Single_Junctions(char* Org_Read,char* Converted_Read,int Read_Skip,int StringLength, unsigned Anchor,int & Inspected_Pairs,SARange & R,int Junc_Count,int Mis_In_Junc_Count,int & Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int Skip,int & Err,int Sign)
{
	assert(StringLength<=READLEN && Converted_Read && Read_Skip>=0 && StringLength >0);
	char Read[200];
	char Read_Head[200];
	PAIR* Pairs= new PAIR[MAX_HITS_TO_STORE+10];if(!Pairs) {cout << "Seek_Junc(): Error allocatind memory..\n";exit(100);}
	Junction Final_Juncs[MAX_JUNCS_TO_STORE+1];Final_Juncs[0].p=UINT_MAX;
	int Label=0;
	LEN L;

	int Fragment=3*MINX;
	if(Read_Skip+2*MINX >StringLength)//cannot penguin 18-mers
	{
		Fragment=MINX+StringLength-Read_Skip;
	}
	B2C(Converted_Read,Read,Fragment);B2C(Org_Read,Read_Head,Read_Skip-MINX);
	L.IGNOREHEAD=0;Split_Read(MINX,L);//we are scanning MINX-mers...
	Inspected_Pairs+=Find_Single_Junc(Read,Read_Head,Converted_Read,Pre,Suf,Fragment,L,Pairs,Final_Juncs,Err,Label,'+',Anchor,R,Skip);
	if (Inspected_Pairs < MAX_INSPECTED_PAIRS) // Dont want to go to an almost infinite loop..
	{
		for(int i=0;Final_Juncs[i].p!=UINT_MAX;i++)
		{
			unsigned p=Final_Juncs[i].p;unsigned q=Final_Juncs[i].q;int r=Final_Juncs[i].r;//assert(q>=p && r>=0);
			//printf("q-p: %d, fragment-2*minx: %d\n", int(q-p), Fragment-2*MINX);
			if(int(q-p)<=(Fragment-2*MINX) || q<p) continue;
			bool Junc_Already_Found=false;
			for(int j=0;j<Compiled_Junctions_Ptr;j++) 
			{
				if(Compiled_Junctions[j].p==Final_Juncs[i].p && Compiled_Junctions[j].q==Final_Juncs[i].q && Compiled_Junctions[j].r==Final_Juncs[i].r )
				{
					Junc_Already_Found=true;break;
				}
			}
			if(Junc_Already_Found) continue;
			Trans_Array[Trans_Array_Ptr]=Final_Juncs[i];//.x=p;Trans_Array[Trans_Array_Ptr].y=q;
			Seek_Junc(Converted_Read+r-Skip,R,0,Junc_Count+1,Mis_In_Junc_Count,q+1,StringLength-r,Trans_Array_Ptr+1,Pre,Suf,Inspected_Pairs,Err,Sign);
			if (Inspected_Pairs >= MAX_INSPECTED_PAIRS) {Err++;break;}
		}
	}
	else Err++;
	delete [] Pairs;
}

int Seek_Single_Strand(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,int Sign)
{
	int c;
	if(Generic_Hits.Lookupsize==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}

	SARange SA;
	SA.Start=0;SA.End=revfmi->textLength;
	SA.Level=1; SA.Skip=0;SA.Mismatches=0;Generic_Hits.Hit_Array_Ptr=0;Generic_Hits.Current_Tag=Current_Tag;

	int Inspected_Pairs=0;
	int Err=0;
	Seek_Junc(Current_Tag,SA,0,0,0,UINT_MAX,StringLength,0,MF_Pre,MF_Suf,Inspected_Pairs,Err,Sign);
	return Err;
}

int Seek_All_Junc(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf)
{
	int Err= Seek_Single_Strand(Current_Tag,StringLength,MF_Pre,MF_Suf,true/*Plus*/);
	char RC_Read[MAXDES],RC_Bin[MAXDES];
	Convert_Reverse(Current_Tag,RC_Read,RC_Bin,StringLength);
	Err+= Seek_Single_Strand(RC_Bin,StringLength,MF_Pre,MF_Suf,false/*Minus*/);
	return Err;
}

void fillProbArray(float prob[16][64][2], string filename) {
	ifstream in;
	in.open(filename.c_str());
	for(int i=0;i<16;i++) {
		for(int j=0;j<64;j++) {
			in >> prob[i][j][0];
			in >> prob[i][j][1];
		}
	}
	in.close();
}

void Init_Prob() {
	fillProbArray(Donor_Prob,"donors.txt");
	fillProbArray(Acc_Prob,"acceptors.txt");
}
