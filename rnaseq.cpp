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
#include <string>
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
#include <pthread.h>
#include "sched.h" 
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
const int MAX_JUNCS_ALLOWED=100;
const int MAX_JUNCS_IN_TRANSCRIPT=20;//Maximum number of junctions that are allowed in a trascript..
const int MAX_INSPECTED_PAIRS=20000;//INT_MAX;

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
unsigned Offsets[80];
Offset_Record Genome_Offsets_Main[80];
int ORG_STRINGLENGTH;

unsigned RQ_Hits;
SA* SA_Index;
int* SA_Blocks;
char COMPRESS;

float Donor_Prob[16][64][2];
float Acc_Prob[16][64][2];
Parameters CL;	
Index_Info Genome_Files;
MMPool *mmPool;
RANGEINDEX Range_Index;
unsigned SOURCELENGTH;
gzFile Input_File,Mate_File;
FILETYPE File_Info;

int THREAD = 0;
pthread_mutex_t lock;
pthread_mutex_t fillProbArraylock;
pthread_mutex_t OpenGenomeFileslock;
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------

struct Transcript_Data
{
	Junction Trans_Array[MAX_JUNCS_IN_TRANSCRIPT+1];
	Junction Compiled_Junctions[MAX_JUNCS_ALLOWED+1];
	Junction Partial_Junctions[MAX_JUNCS_ALLOWED+1];
	int Compiled_Junctions_Ptr;
	int Partial_Junctions_Ptr;
	int Max_Junc_Found;
	int Transcript_Number;
	int Max_Junc_Count;
	int Least_Mis_In_Junc;
	MEMX Generic_Hits;
};

//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*
int Find_All_Single_Junc(char* Read,char* Converted_Read,MEMX & Pre,MEMX & Suf,int STRINGLENGTH,LEN & L,PAIR* Pairs,Junction *Final_Juncs,int & Err,int & Label);
bool Align(char* Source,int StringLength,unsigned Dest_Loc,SARange & R,int Actual_Length,int Read_Skip,Transcript_Data & TD);
int Seek_Junc(char* S,SARange R,int Read_Skip,int Junc_Count,int Mis_In_Junc_Count,unsigned Last_Exon,int StringLength,int Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int & Inspected_Pairs,int & Err,bool Sign,Transcript_Data & TD);
int Seek_Single_Strand(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,int Sign,Transcript_Data & TD);
int Seek_All_Junc(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,Transcript_Data & TD);
int Find_Single_Junc(char* Read,char* Read_Head,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err,int & Label,char Sign,unsigned Anchor,SARange & R,int Skip);
void Enum_Single_Junctions(char* Org_Read,char* Converted_Read,int Read_Skip,int StringLength, unsigned Anchor,int & Inspected_Pairs,SARange & R,int Junc_Count,int Mis_In_Junc_Count,int & Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int Skip,int & Err,int Sign,int Fragment,Transcript_Data & TD);
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
void Reset_MEMX(MEMX & M);
bool Exons_Too_Small(int Trans_Array_Ptr,Transcript_Data & TD);
void *Map(void *T);
void Join_Tables(Offset_Record *Genome_Offsets,int Thread_ID);
void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T);
void Set_Affinity();
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time,Maptime;
	time(&Start_Time);
	if ((pthread_mutex_init(&lock, NULL) != 0) ||(pthread_mutex_init(&fillProbArraylock, NULL) != 0)||(pthread_mutex_init(&OpenGenomeFileslock, NULL) != 0))
	{
		printf("\n mutex init failed\n");
		return 1;
	}

	Parse_Command_line(argc,argv,Genome_Files,CL);
	Load_All_Indexes(Genome_Files,fwfmi,revfmi,mmPool,Range_Index);
	Print_SAM_Header(Annotations,argc,argv,CL.PATTERNFILE);
	Init(revfmi,SOURCELENGTH,Input_File,Mate_File,File_Info,CL,Genome_Files);
	Open_Genome_Files(Genome_Files.LOCATIONFILE,Genome_Offsets_Main,Offsets);
	Conversion_Factor=revfmi->textLength-RQFACTOR;

	Thread_Arg T;
	if (THREAD)
	{
		Launch_Threads(THREAD, Map,T);
	}
	else
	{
		Set_Affinity();
		Map(NULL);
	}
	//Map(NULL);

	UnLoad_Indexes(fwfmi,revfmi,mmPool,Range_Index);
	fprintf(stderr,"\r[++++++++100%%+++++++++]\n");//progress bar....
	time(&End_Time);
	Maptime=difftime(End_Time,Start_Time);

	char Junction_File[]="junctions.bed";
	Print_Junctions(Junction_File,Genome_Offsets_Main);
	time(&End_Time);fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",(float)difftime(End_Time,Start_Time));
	fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",float(Maptime));
	return 0;
}


void *Map(void *T)
{
	Thread_Arg *TA;
	TA=(Thread_Arg*) T;
	int Thread_ID;
	if(THREAD==0 || THREAD==1) {Thread_ID=1;}
	if(T)
	{
		Thread_ID=TA->ThreadID;
	}
	char Str_Thread_ID[20];sprintf(Str_Thread_ID,".%d",Thread_ID);
	Transcript_Data TD;
	unsigned Total_Hits=0,Tags_Processed=0,Tag_Count=0;
	unsigned Number_of_Tags=1000;
	unsigned Progress=0;
	unsigned Hit_ID=0;
	FILE* OUT;
	Offset_Record Genome_Offsets[80];

//------------------- INIT -------------------------------------------------------------
	ofstream SAM[4];
	string SAMFiles[] = {"unique_signal","unique_nosignal",/*"mishits.fq",*/"nonunique_onesignal","others"};
	PAIR *Pairs;
	//Junction Final_Juncs[2*MAX_JUNCS_TO_STORE+2];

	if (!(Pairs=(PAIR*)malloc(sizeof(PAIR)*(MAX_HITS_TO_STORE+10)))) {fprintf(stderr,"Allocate_Memory():malloc error...\n");exit(100);}
	if (CL.JUNCTIONFILE) OUT=File_Open(CL.JUNCTIONFILE,"w"); else OUT=stdout;
	Open_Genome_Files(Genome_Files.LOCATIONFILE,Genome_Offsets,Offsets);
	for(int i=0;i<4;i++)
	{
		Open_Outputs(SAM[i],SAMFiles[i]+Str_Thread_ID+".sam");
	}
	ofstream rejectedSAM;
	string rejectedSAMFile = "rejected";
	Open_Outputs(rejectedSAM, rejectedSAMFile+Str_Thread_ID+".sam");
	
//------------------- INIT -------------------------------------------------------------

	READ Head,Tail;
	int LOOKUPSIZE=3;
	MEMLOOK MLook;MLook.Lookupsize=3;
	Build_Tables(fwfmi,revfmi,MLook);
	LEN L;L.IGNOREHEAD=0;
	Split_Read(RQFACTOR,L);//we are scanning 18-mers...
//--------------------- Setup Data Structure for Batman Prefix ----------------------------------------
	MEMX MF_Pre,MF_Suf;//MemX is the data structure for doing Batman alignment. MF_Pre is for the prefix, MC for suffix..
	Init_Batman(MF_Pre,L,MLook,MAX_MISMATCHES);
	Init_Batman(MF_Suf,L,MLook,MAX_MISMATCHES);
	Init_Batman(TD.Generic_Hits,L,MLook,MAX_MISMATCHES);
//--------------------- Setup Data Structure for Batman End----------------------------------------
	

//--------------------- Load probability information --------------------------
	Init_Prob();

	fprintf(stderr,"======================]\r[");//progress bar....
	int Actual_Tag=0;
	int selectedJunctions[MAX_JUNCS_ALLOWED+1];
	while (Read_Tag(Head,Tail,Input_File,Mate_File,File_Info))
	{
		if(Thread_ID==1 && !Progress_Bar(CL,Number_of_Tags,Progress,Tag_Count,File_Info)) break;
		if(CL.MAX_TAGS_TO_PROCESS && CL.MAX_TAGS_TO_PROCESS<Actual_Tag) break;

		int Label=0;
		Actual_Tag++;

		TD.Max_Junc_Found=INT_MAX;
		TD.Least_Mis_In_Junc=INT_MAX;
		TD.Compiled_Junctions_Ptr=0;
		TD.Partial_Junctions_Ptr=0;
		TD.Transcript_Number=0;
		TD.Max_Junc_Count=0;
		int Err=Seek_All_Junc(Head.Tag,File_Info.STRINGLENGTH,MF_Pre,MF_Suf,TD);
		TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr].p=UINT_MAX;

		if((TD.Compiled_Junctions_Ptr) && !(DEBUG && Err))
		{
			int firstSignal = -2;
			int tempType = Classify_Hits(TD.Compiled_Junctions,firstSignal);
			int approvedPtr;
			approvedPtr = getBest(Head.Tag,File_Info.STRINGLENGTH,TD.Compiled_Junctions, selectedJunctions, true);

			if(approvedPtr > 1) 
			{
				Hit_ID++;
				for(int i=0; i<approvedPtr; i++) 
				{
					Print_Hits(Head,TD.Compiled_Junctions,OUT,rejectedSAM,Tag_Count,selectedJunctions[i],tempType,Hit_ID,Err,Genome_Offsets);//tempType);
				}
			}	
			else 
			{
				//assert(approvedPtr);
				for(int i=0; i<approvedPtr; i++) {
					Print_Hits(Head,TD.Compiled_Junctions,OUT,SAM[tempType],Tag_Count,selectedJunctions[i],tempType,0,Err,Genome_Offsets);//tempType);
				}
			}

		}
	}
	Join_Tables(Genome_Offsets,Thread_ID);

}

void Join_Tables(Offset_Record *Genome_Offsets,int Thread_ID)
{
	pthread_mutex_lock(&lock);
	cout <<"Joining junctions from Thread " << Thread_ID<<endl; 

	for(int i=0;Genome_Offsets[i].Offset !=INT_MAX;i++)
	{
		OP JPair;
		JStat JStat;
		Junction J;

		Hash *Junctions=Genome_Offsets[i].Junc_Hash,*Main=Genome_Offsets_Main[i].Junc_Hash;
		char Junc_Not_Empty = Junctions->Init_Iterate(JPair,JStat);

		while (Junc_Not_Empty)
		{

			J.Type=JStat.Junc_Type;
			Main->Insert(JPair,J,JStat.Unique);
			Junc_Not_Empty=Junctions->Iterate(JPair,JStat);
		}
	}
	cout <<"Done.."<<endl;
	pthread_mutex_unlock(&lock);
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



int Classify_Hits(Junction * Final_Juncs, int & firstSignal){
	int signalCount = 0, count = 0;
	for(int i=0; Final_Juncs[i].p != UINT_MAX; i++){
		count++;
		if(Final_Juncs[i].q)//!=Final_Juncs[i].p-1)//not an exact match?
		{
			Final_Juncs[i].Type=Final_Juncs[i].isCanonical();
			if(Final_Juncs[i].Type) {
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
	int Mismatch_Pos[MAXDES];
	Convert_Reverse(Current_Tag,RC_Read,RC_Bin,StringLength);

//Calculate true mismatches..
	int Junc_Count=Final_Juncs[0].Junc_Count;
	for(int i =0; Final_Juncs[i].p != UINT_MAX;i+=(Junc_Count ? Junc_Count:1),\
						   Junc_Count=Final_Juncs[i].Junc_Count) 
	{
		if(Final_Juncs[i].q)
		{
			int Last=0;
			for(int j=0;j<Junc_Count;j++)
			{
				Get_Bases_ASCII(Final_Juncs[i+j].p-(Final_Juncs[i+j].r), (Final_Juncs[i+j].r), Cat+Last);
				Last+=Final_Juncs[i+j].r;
			}
			Get_Bases_ASCII(Final_Juncs[i+Junc_Count-1].q+1, StringLength-Last, Cat+Last);
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

		int Pos=0;
		for(int j=0;j<StringLength;j++)
		{
			if(Cat[j]!="ACGT"[Read[j]]) 
			{
				Final_Juncs[i].Mismatches++;
				Mismatch_Pos[Pos++]=j;
			}
			Cat[j]=0;
		}

		Final_Juncs[i].score= 0;
		/*if(Final_Juncs[i].r>(StringLength-Final_Juncs[i].r))//large part is in first anchor..
		{
			if(Mismatch_Pos[0]>Final_Juncs[i].r)
			{
				Final_Juncs[i].score= -1;
			}
		}
		else
		{
			if(Mismatch_Pos[0]<Final_Juncs[i].r)
			{
				Final_Juncs[i].score= -1;
			}
		}
		Final_Juncs[i].score= Mismatch_Pos[0];*/
	}

	for(int i =0; Final_Juncs[i].p != UINT_MAX; i+=(Final_Juncs[i].Junc_Count)? Final_Juncs[i].Junc_Count:1)
	{
		int Junc_Score=0;
		if(!Final_Juncs[i].q) Junc_Score=4;
		else
		{
			int Junc_Count=Final_Juncs[i].Junc_Count;
			int Label=Final_Juncs[i].Label;
			unsigned dist= 0;

			for(int j=0;j<Junc_Count;j++)
			{
				assert(Label==Final_Juncs[i+j].Label);
				Junc_Score+=Final_Juncs[i+j].Type;
				if(dist<(Final_Juncs[i+j].q-Final_Juncs[i+j].p))
				{
					dist=Final_Juncs[i+j].q-Final_Juncs[i+j].p;
				}
			}

			assert(dist>0);
			if(dist>30000) Junc_Score+=1;
			else if(dist>5000) Junc_Score+=2;
			else if(dist>2000) Junc_Score+=3;
			else if(dist>=100) Junc_Score+=4;
			else if(dist>60) Junc_Score+=1;
		}
		int tempScore = -3*Final_Juncs[i].Mismatches+Junc_Score;//+Final_Juncs[i].score;
		if(Final_Juncs[i].Junc_Count<=4 && tempScore >= max)
		{
			if(tempScore!=max) ptr = 0;
			max = tempScore;
			approved[ptr] = i;
			++ptr;
		}

	}
	//assert(ptr);
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
		//TODO: enable this for speed
		/*if(Skip==0)// || (R.End-R.Start < SAGAP_CUTOFF_T))//first l-mer..
		{
			Reset_MEMX(MF_Pre);
			MF_Pre.Hit_Array[0]=R;MF_Pre.Hit_Array_Ptr=1;
			MF_Pre.Exact_Match_Forward[L.LH-1].Start=0;//This is workaround.. need to store the exact value..
		}
		else*/
		{
			MF_Pre.Current_Tag=Converted_Read;MF_Pre.Hit_Array_Ptr=0;MF_Pre.Hits=0;
			Last_Mis_Pre=Scan(MF_Pre,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
		}
	}
	else
	{
		MF_Pre.Hit_Array[0].Start=MF_Pre.Hit_Array[0].End=Anchor;MF_Pre.Hit_Array_Ptr=1;
		MF_Pre.Hit_Array[1].Start=0;
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
					assert(Skip>=MINX || Skip==0);
					junctions[i].r+=Skip;
					S.p=junctions[i].p;
					S.q=junctions[i].q; 	
					S.r=junctions[i].r;

					junctions[i].Chrom=A.Name;
					junctions[i].Label=Label;
					junctions[i].Sign=Sign;
					junctions[i].ID=A.ID;

					/*FILE* F=Genome_Offsets[A.ID].Out_File;
					fwrite(&S,sizeof(Split_Map),1,F);*/

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

bool Align(char* Source,int StringLength,unsigned Dest_Loc,SARange & R,int Actual_Length,int Read_Skip,Transcript_Data & TD)
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
			TD.Generic_Hits.Hit_Array[0]=R;
			TD.Generic_Hits.Hit_Array_Ptr=1;
			return true;
		}
	}	
	else
	{
		Reset_MEMX(TD.Generic_Hits);TD.Generic_Hits.Hit_Array_Ptr=0;R.Level=1;TD.Generic_Hits.FMIndex=REVERSE;
		int T_Len=TD.Generic_Hits.L.STRINGLENGTH;TD.Generic_Hits.L.STRINGLENGTH=Actual_Length;
		int Lo=0;
		Extend_Forwards(Source,R,/*MAX_MISMATCHES*/MIS_DENSITY,1,StringLength,INT_MAX,revfmi,TD.Generic_Hits,true,Read_Skip,Lo);
		if(Lo)
		{
			SARange SA;
			SA=TD.Generic_Hits.Hit_Array[Lo];
			TD.Generic_Hits.Hit_Array[Lo]=TD.Generic_Hits.Hit_Array[0];
			TD.Generic_Hits.Hit_Array[0]=SA;
			assert(SA.Mismatches==R.Mismatches);
		}
		TD.Generic_Hits.L.STRINGLENGTH=T_Len;
		if(TD.Generic_Hits.Hit_Array_Ptr) return true;
		else return false;
	}
}

inline void B2C(char *S,char *D,int StringLength)
{
	*(D+StringLength)=0;
	for(int i=0;i<StringLength;i++,S++,D++)
		*D="ACGT"[*S];
}

int Seek_Junc(char* S,SARange R,int Read_Skip,int Junc_Count,int Mis_In_Junc_Count,unsigned Last_Exon,int StringLength,int Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int & Inspected_Pairs,int & Err,bool Sign,Transcript_Data & TD)
{
	assert(Trans_Array_Ptr<MAX_JUNCS_IN_TRANSCRIPT && TD.Compiled_Junctions_Ptr>=0);assert(Read_Skip>=0 && StringLength >=0 && StringLength<=READLEN && Inspected_Pairs>=0);
	assert(Read_Skip>=MINX || Read_Skip==0);
	if(Read_Skip>StringLength-MINX)//End of read..
	{
		if(Align(S+Read_Skip,StringLength-Read_Skip,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip,R,StringLength,Read_Skip,TD))//End of the read can be matched..
		{
			if(!Junc_Count)//Match..
			{
				TD.Max_Junc_Found=Junc_Count;
				Junction J;J.Sign=(Sign)? 1:0;J.Mismatches=R.Mismatches;J.Junc_Count=0;
				if (Last_Exon!=UINT_MAX)//one hit.. 
				{
					J.p=Last_Exon;
					J.r=J.q=0;
					if(TD.Compiled_Junctions_Ptr<MAX_JUNCS_ALLOWED)
					{
						TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr++]=J;
					}
					else 
					{
						Err=1;
					}
				}
				else//Many hits..
				{
					for(unsigned i=R.Start;i<=R.End;i++)
					{
						J.p=revfmi->textLength-READLEN-BWTSaValue(revfmi,i);;
						J.r=J.q=0;
						if(TD.Compiled_Junctions_Ptr<MAX_JUNCS_ALLOWED)
						{
							TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr++]=J;
						}
						else 
						{
							Err=1;break;
						}
					}
				}
				assert(TD.Compiled_Junctions_Ptr>0 && TD.Compiled_Junctions_Ptr<=MAX_JUNCS_ALLOWED);
				return TRANCRIPT_END;
			}

			if(Junc_Count) 
			{
				TD.Max_Junc_Found=Junc_Count;
			}
			if(TD.Least_Mis_In_Junc> Mis_In_Junc_Count) TD.Least_Mis_In_Junc=Mis_In_Junc_Count; 
			TD.Transcript_Number++;
			if(TD.Max_Junc_Count<Junc_Count)//count max transcribed juncs..
				TD.Max_Junc_Count=Junc_Count;

			if(TD.Compiled_Junctions_Ptr+Junc_Count<= MAX_JUNCS_ALLOWED)
			{
				for(int i=0;i<Junc_Count;i++)
				{
					//assert(Junc_Count<=2);
					TD.Trans_Array[i].Junc_Count=Junc_Count;TD.Trans_Array[i].Label=TD.Transcript_Number;TD.Trans_Array[i].Sign=(Sign)? 1:0;TD.Trans_Array[i].Mismatches=R.Mismatches;TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr++]=TD.Trans_Array[i];
				}
			}
			else
			{
				Err=1;
			}

			assert(TD.Compiled_Junctions_Ptr>0 && TD.Compiled_Junctions_Ptr<=MAX_JUNCS_ALLOWED);
			return TRANCRIPT_END;
		}
		else//End of the read may contain extra junc that cannot be processed..
		{
			if(Junc_Count && TD.Least_Mis_In_Junc> Mis_In_Junc_Count) TD.Least_Mis_In_Junc=Mis_In_Junc_Count; 

			int JC=(Junc_Count)? Junc_Count:1;//TODO
			if(TD.Partial_Junctions_Ptr+JC< MAX_JUNCS_ALLOWED)
			{
				for(int i=0;i<JC;i++)
				{
					TD.Trans_Array[i].Junc_Count=Junc_Count;TD.Trans_Array[i].Label=TD.Transcript_Number;TD.Trans_Array[i].Sign=(Sign)? 1:0;TD.Trans_Array[i].Mismatches=R.Mismatches;TD.Partial_Junctions[TD.Partial_Junctions_Ptr++]=TD.Trans_Array[i];
				}
			}
			else
			{
				Err=1;
			}

			assert(TD.Partial_Junctions_Ptr>0 && TD.Partial_Junctions_Ptr<=MAX_JUNCS_ALLOWED);
			return DUMMY_JUNC;
		}
	}

	if(Align(S+Read_Skip,MINX,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip,R,Read_Skip+MINX,Read_Skip,TD))//Can exact extension be done..
	{
		std::vector <SARange> Hits(TD.Generic_Hits.Hit_Array_Ptr);
		for (int i=0;i<TD.Generic_Hits.Hit_Array_Ptr;i++)//Save hits in a vector..
		{
			assert(TD.Generic_Hits.Hit_Array[i].Start);
			Hits[i]=TD.Generic_Hits.Hit_Array[i];
		}
		unsigned Last_ExonT;int Hit_Count=TD.Generic_Hits.Hit_Array_Ptr-1;
		for (int i=0;i<=Hit_Count;i++)//Save hits in a vector..
		{
			SARange SA=Hits[i];
			Last_ExonT=Last_Exon;
			if(SA.Start==SA.End && Last_Exon==UINT_MAX) 
			{
				Last_ExonT=SA.Start-Read_Skip;
			}
			Seek_Junc(S,SA,Read_Skip+MINX,Junc_Count,Mis_In_Junc_Count,Last_ExonT,StringLength,Trans_Array_Ptr,Pre,Suf,Inspected_Pairs,Err,Sign,TD);
			if(Err) return 0;
			if (Inspected_Pairs >= MAX_INSPECTED_PAIRS) {Err++;return 0;}
		}
	}

//Junction check routines
	if(Read_Skip>=MINX)
	{
		if(Exons_Too_Small(Trans_Array_Ptr,TD)) return 0;
		//Case of a single junction between S[a,a+2l-1]
		if (Junc_Count-1<TD.Max_Junc_Found) //This will add one junction..
		{
			//the junction is between S[Read_Skip,Read_Skip+18]
			//so anchor S[Read_Skip-18,Read_Skip]
			//Penguin S[Read_Skip-18,Read_Skip] ---- S[Read_Skip+18,Read_Skip+36] and extend..
			Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX,TD);
			//Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX-6,TD);
			//Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX-12,TD);
			//Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX-15,TD);
		}
	}

	return Trans_Array_Ptr;

}

void Enum_Single_Junctions(char* Org_Read,char* Converted_Read,int Read_Skip,int StringLength, unsigned Anchor,int & Inspected_Pairs,SARange & R,int Junc_Count,int Mis_In_Junc_Count,int & Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int Skip,int & Err,int Sign,int Fragment,Transcript_Data & TD)
{
	assert(StringLength<=READLEN && Converted_Read && Read_Skip>=0 && StringLength >0);
	char Read[200];
	char Read_Head[200];
	PAIR* Pairs= new PAIR[MAX_HITS_TO_STORE+10];if(!Pairs) {cout << "Seek_Junc(): Error allocatind memory..\n";exit(100);}
	Junction Final_Juncs[MAX_JUNCS_TO_STORE+1];Final_Juncs[0].p=UINT_MAX;
	int Label=0;
	LEN L;

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
			for(int j=0;j<TD.Compiled_Junctions_Ptr;j++) 
			{
				if(TD.Compiled_Junctions[j].p==Final_Juncs[i].p && TD.Compiled_Junctions[j].q==Final_Juncs[i].q && TD.Compiled_Junctions[j].r==Final_Juncs[i].r )
				{
					Junc_Already_Found=true;break;
				}
			}
			if(Junc_Already_Found) continue;
			TD.Trans_Array[Trans_Array_Ptr]=Final_Juncs[i];//.x=p;Trans_Array[Trans_Array_Ptr].y=q;
			Seek_Junc(Converted_Read+r-Skip,R,0,Junc_Count+1,Mis_In_Junc_Count,q+1,StringLength-r,Trans_Array_Ptr+1,Pre,Suf,Inspected_Pairs,Err,Sign,TD);
			if(Err) break;
			if (Inspected_Pairs >= MAX_INSPECTED_PAIRS) {Err++;break;}
		}
	}
	else Err++;
	delete [] Pairs;
}

int Seek_Single_Strand(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,int Sign,Transcript_Data & TD)
{
	SARange SA;
	SA.Start=0;SA.End=revfmi->textLength;//SA.FMIndex=REVERSE;
	SA.Level=1; SA.Skip=0;SA.Mismatches=0;SA.Mismatch_Char=0;Reset_MEMX(TD.Generic_Hits);TD.Generic_Hits.Hit_Array_Ptr=0;TD.Generic_Hits.Current_Tag=Current_Tag;

	int Inspected_Pairs=0;
	int Err=0;
	Seek_Junc(Current_Tag,SA,0,0,0,UINT_MAX,StringLength,0,MF_Pre,MF_Suf,Inspected_Pairs,Err,Sign,TD);
	return Err;
}

int Seek_All_Junc(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,Transcript_Data & TD)
{
	int Err= Seek_Single_Strand(Current_Tag,StringLength,MF_Pre,MF_Suf,1/*Plus*/,TD);
	if (Err) return Err;
	char RC_Read[MAXDES],RC_Bin[MAXDES];
	Convert_Reverse(Current_Tag,RC_Read,RC_Bin,StringLength);
	Err+= Seek_Single_Strand(RC_Bin,StringLength,MF_Pre,MF_Suf,0/*Minus*/,TD);
	return Err;
}

void fillProbArray(float prob[16][64][2], string filename) {
	ifstream in;
	pthread_mutex_lock(&fillProbArraylock);
	in.open(filename.c_str());
	for(int i=0;i<16;i++) {
		for(int j=0;j<64;j++) {
			in >> prob[i][j][0];
			in >> prob[i][j][1];
		}
	}
	in.close();
	pthread_mutex_unlock(&fillProbArraylock);
}

void Init_Prob() {
	fillProbArray(Donor_Prob,"donors.txt");
	fillProbArray(Acc_Prob,"acceptors.txt");
}

void Reset_MEMX(MEMX & M)
{
	M.Left_Mishits_Pointer=0;
	M.Right_Mishits_Pointer=0;
	M.Possible_20_Pointer=0;
	M.Possible_02_Pointer=0;
	M.Mismatches_Forward_Pointer=0;
	M.Mismatches_Backward_Pointer=0;
	M.Two_Mismatches_At_End_Pointer=0;
	M.Two_Mismatches_At_End_Forward_Pointer=0;
	M.Hit_Array_Ptr=0;
	M.Hits=0;
}

bool Exons_Too_Small(int Trans_Array_Ptr,Transcript_Data & TD)
{
	int Count=0;
	for(int i=0;i<Trans_Array_Ptr;i++)
	{
		if(TD.Trans_Array[i].r<=20) Count++;
		if(TD.Trans_Array[i].r<=18) return true;
	}
	if(Count>1)
		return true;
	else
		return false;
}


void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T)
{
	Threading* Thread_Info=(Threading*) malloc(sizeof(Threading)*NTHREAD);
	int Thread_Num=0;
	pthread_attr_t Attrib;
	pthread_attr_init(&Attrib);
	pthread_attr_setdetachstate(&Attrib, PTHREAD_CREATE_JOINABLE);

	for (int i=0;i<NTHREAD;i++)
	{
		T.ThreadID=i;
		Thread_Info[i].Arg=T;
		//if(!(Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,NULL,Map_t,(void*) &Thread_Info[i].Arg))) Thread_Num++;
		Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,&Attrib,Map_t,(void*) &Thread_Info[i].Arg);
		if(Thread_Info[i].r) {printf("Launch_Threads():Cannot create thread..\n");exit(-1);} else Thread_Num++;
	}
	printf("%d Threads runnning ...\n",Thread_Num);
	pthread_attr_destroy(&Attrib);

	for (int i=0;i<NTHREAD;i++)
	{
		pthread_join(Thread_Info[i].Thread,NULL);
	}
}

void Set_Affinity()
{
	cpu_set_t Set;
	CPU_ZERO(&Set);

	if(sched_getaffinity(0,sizeof(cpu_set_t),&Set)<0)
	{
		printf("Affinity could not be get..\n");
	}
	else
	{
		for (int i=0;i<CPU_SETSIZE;i++)
		{
			if(CPU_ISSET(i,&Set))
			{
				printf("Bound to %d\n",i);
				CPU_ZERO(&Set);
				CPU_SET(i, &Set);
				if(sched_setaffinity(0, sizeof(Set), &Set)<0)
				{
					printf("Affinity could not be set..\n");
				}
				return;
			}
		}
	}

}
