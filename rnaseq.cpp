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
#define QLIMIT_FLOAT 30.0f
#define QLIMIT QLIMIT_FLOAT 

extern const int UNIQUE_SIGNAL=0;
extern const int UNIQUE_NOSIGNAL=1;
extern const int NON_UNIQUE_SIGNAL=2;
extern const int NON_UNIQUE_OTHERS=3;
//extern const int MINX=18;//Minimum extension..
extern const int EXON_GEN_LEN=5000;
extern const int CANONICAL_SCORE=3;

const int DUMMY_JUNC=1000;
const int TRANCRIPT_END=1001;
extern const bool DEBUG=false;
int QSUM_LIMIT=32;
//const int QSUM_LIMIT=60;

struct Transcript
{
	unsigned x;
	unsigned y;
};

//Transcript Trans_Array[20];
const int MAX_JUNCS_ALLOWED=100;//00;
const int MAX_JUNCS_IN_TRANSCRIPT=20;//Maximum number of junctions that are allowed in a trascript..
const int MAX_INSPECTED_PAIRS=20000;//INT_MAX;
const int QUALITYCONVERSIONFACTOR=33;
bool ONEMULTIHIT=true;
bool SOFTCLIP=true;

const int POWLIMIT=300;
bool DEEPSCAN=false;
char Dummy_Q[]="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
float POW10[POWLIMIT];
bool PRINT_NON_CANON=false;
bool PRE_MAP=true;

int TENMER=10;
int MIS_DENSITY=1;
int RQFACTOR=18;
int MINX=18;
int INDEX_RESOLUTION=30000;
int EXONGAP;
int READLEN;
bool SAM_READER=true;
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
int INIT_MIS_SCAN=-1;

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
FILE *Input_File,*Mate_File;
FILETYPE File_Info;

int THREAD = 0;
pthread_mutex_t lock;
pthread_mutex_t fillProbArraylock;
pthread_mutex_t OpenGenomeFileslock;
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------
struct Potential_Hit
{
	SARange SA;
	int Read_Skip;
	char Sign;
};

struct Transcript_Data
{
	Junction Trans_Array[MAX_JUNCS_IN_TRANSCRIPT+1];
	Junction Compiled_Junctions[MAX_JUNCS_ALLOWED+1];
	Junction Partial_Junctions[MAX_JUNCS_ALLOWED+1];
	Potential_Hit Ext_Array[MAX_JUNCS_ALLOWED+1];
	int Compiled_Junctions_Ptr;
	int Partial_Junctions_Ptr;
	int Max_Junc_Found;
	int Transcript_Number;
	int Ext_Array_Ptr;
	int Max_Junc_Count;
	int Least_Mis_In_Junc;
	MEMX Generic_Hits;
};

//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*
int Find_All_Single_Junc(char* Read,char* Converted_Read,MEMX & Pre,MEMX & Suf,int STRINGLENGTH,LEN & L,PAIR* Pairs,Junction *Final_Juncs,int & Err,int & Label);
bool Align(char* Source,int StringLength,unsigned Dest_Loc,SARange & R,int Actual_Length,int Read_Skip,Transcript_Data & TD,char* Q);
int Seek_Junc(char* S,SARange R,int Read_Skip,int Junc_Count,int Mis_In_Junc_Count,unsigned Last_Exon,int StringLength,int Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int & Inspected_Pairs,int & Err,bool Sign,Transcript_Data & TD,char* Q);
int Seek_Single_Strand(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,int Sign,Transcript_Data & TD,char *Q);
int Seek_All_Junc(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,Transcript_Data & TD,char* Q);
int Find_Single_Junc(char* Read,char* Read_Head,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err,int & Label,char Sign,unsigned Anchor,SARange & R,int Skip);
void Enum_Single_Junctions(char* Org_Read,char* Converted_Read,int Read_Skip,int StringLength, unsigned Anchor,int & Inspected_Pairs,SARange & R,int Junc_Count,int Mis_In_Junc_Count,int & Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int Skip,int & Err,int Sign,int Fragment,Transcript_Data & TD,char* Q);
void Load_All_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,RANGEINDEX & Range_Index);
void Init(BWT *revfmi,unsigned & SOURCELENGTH,PAIR* & Pairs,FILE* & Input_File,FILE* & Mate_File,FILETYPE & File_Info,Parameters & CL,FILE* & OUT,Index_Info & Genome_Files);
bool  Progress_Bar(Parameters & CL,unsigned & Number_of_Tags,unsigned & Progress,unsigned & Tag_Count,FILETYPE & File_Info);
int Find_Pairings(int & Pairs_Index,SARange* MF_Pre_Hits,SARange* MF_Suf_Hits,PAIR* &  Pairs, char* Read,char* Read_Head,int STRINGLENGTH,int & Final_Juncs_Ptr,Junction *Final_Juncs,int & Min_Mismatch,int Mis_In_Head,int Mis_In_Tail,int Skip,int & Label,char Sign);
void Open_Outputs(ofstream & SAM,string filename);
inline char* Nullify_String(char* S);
//void Print_Hits(READ & Head,Junction *Final_Juncs,FILE* OUT,ofstream & SAM,int Tag_Count, int firstSignal);
int Classify_Hits(Junction * Final_Juncs, int & firstSignal);
inline void Convert_Reverse(char* Read,char * RC_Read,char* RC_Bin,int StringLength);
void Print_SAM_Header(std::map <unsigned, Ann_Info> Annotations,int argc,char* argv[],char* Input_File);
void Init_Prob();
void fillProbArray(float prob[16][64][2], string filename);
//int getBest(char* Current_Tag,int StringLength,Junction * Compiled_Junctions, int * approved, bool print);
int getBest(char* Current_Tag,int StringLength,Junction * Final_Juncs, int * approved, bool print,char* Q);
void Reset_MEMX(MEMX & M);
bool Exons_Too_Small(int Trans_Array_Ptr,Transcript_Data & TD);
void *Map(void *T);
bool Messy_CIGAR(char *Cig,int StringLength);
void Join_Tables(Offset_Record *Genome_Offsets,int Thread_ID);
void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T);
void Set_Affinity();
void Scan_End_Junc(char* Current_Tag,int StringLength,Transcript_Data & TD,int & Err,char* Q);
void Do_End_Scan(char* Current_Tag,char* Current_Tag_ASCII,int StringLength,char Sign,Transcript_Data & TD,int & Err,char* Q);
float Penalty(char Q);
void Build_Pow10();
void Reverse_Quality(char* Dest,char* Q,int StringLength);
bool Mismatch_Hit_Nice(int Mismatch_Scan,MEMX & MF,MEMX & MC,char* Q,int StringLength);
void Scan_Partial_Read(char* Current_Tag,int StringLength,MEMX & Pre, MEMX & Suf,Transcript_Data & TD,char* Q);
void Partial_Scan(char* Current_Tag,char* Current_Tag_ASCII,int StringLength,MEMX & Pre,MEMX & Suf,char Sign,Transcript_Data & TD,int & Err,char* Q);
void Scan_Right_Read(char* Current_Tag,int StringLength,MEMX & Pre, MEMX & Suf,Transcript_Data & TD,char* Q);
void Extension_Scan(char* Current_Tag,char* Current_Tag_ASCII,int StringLength,MEMX & Pre,MEMX & Suf,char Sign,Transcript_Data & TD,int & Err,char* Q);
void Seek_Residue(int & GT_Ptr,int *GT,char* Current_Tag_ASCII,int StringLength,char* RightX,char* Q,char *Signal,unsigned Loc,Transcript_Data & TD,char Sign);
void Print_Matches(MEMX & MF,MEMX & MC,READ & Head,ofstream & SAM,int StringLength,Offset_Record *Genome_Offsets,int Mismatches);
void Print_Unmapped(READ & Head,int StringLength,ofstream & MISHIT);
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

	Build_Pow10();
	Parse_Command_line(argc,argv,Genome_Files,CL);
	Load_All_Indexes(Genome_Files,fwfmi,revfmi,mmPool,Range_Index);
	Print_SAM_Header(Annotations,argc,argv,CL.PATTERNFILE);
	Init(revfmi,SOURCELENGTH,Input_File,Mate_File,File_Info,CL,Genome_Files,INIT_MIS_SCAN);
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
	ofstream SAM[5];
	string SAMFiles[] = {"alignments","unique_nosignal",/*"mishits.fq",*/"nonunique_onesignal","others","mishits"};
	PAIR *Pairs;
	//Junction Final_Juncs[2*MAX_JUNCS_TO_STORE+2];

	if (!(Pairs=(PAIR*)malloc(sizeof(PAIR)*(MAX_HITS_TO_STORE+10)))) {fprintf(stderr,"Allocate_Memory():malloc error...\n");exit(100);}
	if (CL.JUNCTIONFILE) OUT=File_Open(CL.JUNCTIONFILE,"w"); else OUT=stdout;
	Open_Genome_Files(Genome_Files.LOCATIONFILE,Genome_Offsets,Offsets);
	//for(int i=0;i<5;i++)
	for(int i=0;i<1;i++)
	{
		Open_Outputs(SAM[i],SAMFiles[i]+Str_Thread_ID+".sam");
	}
	/*ofstream rejectedSAM;
	string rejectedSAMFile = "rejected";
	Open_Outputs(rejectedSAM, rejectedSAMFile+Str_Thread_ID+".sam");*/
	
//------------------- INIT -------------------------------------------------------------

	READ Head,Tail;
	SAMREAD SAM_Read;SAM_Read.Enable=SAM_READER;
	int LOOKUPSIZE=3;
	MEMLOOK MLook;MLook.Lookupsize=3;
	MEMLOOK MLook_Long;MLook_Long.Lookupsize=6;
	Build_Tables(fwfmi,revfmi,MLook);
	Build_Tables(fwfmi,revfmi,MLook_Long);
	LEN L,L_Long;L.IGNOREHEAD=L_Long.IGNOREHEAD=0;
	Split_Read(RQFACTOR,L);//we are scanning 18-mers...
	Split_Read(File_Info.STRINGLENGTH,L_Long);//we are scanning whole read...
//--------------------- Setup Data Structure for Batman Prefix ----------------------------------------
	MEMX MF_Pre,MF_Suf,MF,MC;//MemX is the data structure for doing Batman alignment. MF_Pre is for the prefix, MC for suffix..
	Init_Batman(MF_Pre,L,MLook,MAX_MISMATCHES);
	Init_Batman(MF_Suf,L,MLook,MAX_MISMATCHES);
	Init_Batman(TD.Generic_Hits,L,MLook,MAX_MISMATCHES+1);

	Init_Batman(MF,L_Long,MLook_Long,MAX_MISMATCHES);
	Init_Batman(MC,L_Long,MLook_Long,MAX_MISMATCHES);
//--------------------- Setup Data Structure for Batman End----------------------------------------
	

//--------------------- Load probability information --------------------------
	Init_Prob();

	fprintf(stderr,"======================]\r[");//progress bar....
	int Actual_Tag=0;
	int selectedJunctions[MAX_JUNCS_ALLOWED+1];
	if(SAM_READER) PRE_MAP=false;
	while (Read_Tag(Head,Tail,Input_File,Mate_File,File_Info,SAM_Read))
	{
		if(Thread_ID==1 && !Progress_Bar(CL,Number_of_Tags,Progress,Tag_Count,File_Info)) break;
		if(CL.MAX_TAGS_TO_PROCESS && CL.MAX_TAGS_TO_PROCESS<Actual_Tag) break;

		int Label=0;
		Actual_Tag++;
		//cout<<"-------------------"<<Actual_Tag<<"--------------------------"<<endl;

		TD.Max_Junc_Found=INT_MAX;
		TD.Least_Mis_In_Junc=INT_MAX;
		TD.Compiled_Junctions_Ptr=0;
		TD.Partial_Junctions_Ptr=0;
		TD.Ext_Array_Ptr=0;
		TD.Transcript_Number=0;
		TD.Max_Junc_Count=0;


		char Rev_Bin[MAXTAG],Rev[MAXTAG];
		int READLEN_T=File_Info.STRINGLENGTH,Skip_Bases=0,Real_Skip=0;
		Convert_Reverse(Head.Tag,Rev,Rev_Bin,File_Info.STRINGLENGTH);
		MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Current_Tag=Head.Tag;MF_Pre.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
		MC.Hits=0;MC.Hit_Array_Ptr=0;MC.Current_Tag=Rev_Bin;MC.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
		int Mismatch_Scan= -1;//Scan_Both(MF,MC,INIT_MIS_SCAN,L_Long,fwfmi,revfmi,0,2);
		/*if(PRE_MAP && Mismatch_Hit_Nice(Mismatch_Scan,MF,MC,(File_Info.FILETYPE==FQ)? Head.Quality:Dummy_Q,File_Info.STRINGLENGTH)) 
		{
			Print_Matches(MF,MC,Head,SAM[0],READLEN_T,Genome_Offsets,Mismatch_Scan);
			continue;
		}*/
		if(SAM_Read.NM<1 && !Messy_CIGAR(SAM_Read.Cigar,File_Info.STRINGLENGTH))
		{
			SAM[0] << SAM_Read.SAM_Line;
			continue;
		}

		Head.Tag_Copy[READLEN_T]=0;
		int Err=Seek_All_Junc(Head.Tag,READLEN_T,MF_Pre,MF_Suf,TD,(File_Info.FILETYPE==FQ)? Head.Quality:Dummy_Q);
		TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr].p=UINT_MAX;

		if(!TD.Compiled_Junctions_Ptr)//try junc in the end..
		{
			Scan_End_Junc(Head.Tag,READLEN_T,TD,Err,(File_Info.FILETYPE==FQ)? Head.Quality:Dummy_Q);

			if(!TD.Compiled_Junctions_Ptr)
			{
				for(int i=1;i<File_Info.STRINGLENGTH-50;i+=5)
				{
					TD.Max_Junc_Found=INT_MAX;
					TD.Least_Mis_In_Junc=INT_MAX;
					TD.Compiled_Junctions_Ptr=0;
					TD.Partial_Junctions_Ptr=0;
					TD.Ext_Array_Ptr=0;
					TD.Transcript_Number=0;
					TD.Max_Junc_Count=0;

					READLEN_T-=5;Skip_Bases+=5;
					Err=Seek_All_Junc(Head.Tag+Skip_Bases,READLEN_T,MF_Pre,MF_Suf,TD,(File_Info.FILETYPE==FQ)? Head.Quality+Skip_Bases:Dummy_Q);
					TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr].p=UINT_MAX;
					if(TD.Compiled_Junctions_Ptr)
					{
						Real_Skip=Skip_Bases;
						break;
					}
					TD.Max_Junc_Found=INT_MAX;
					TD.Least_Mis_In_Junc=INT_MAX;
					TD.Compiled_Junctions_Ptr=0;
					TD.Partial_Junctions_Ptr=0;
					TD.Ext_Array_Ptr=0;
					TD.Transcript_Number=0;
					TD.Max_Junc_Count=0;
					Err=Seek_All_Junc(Head.Tag,READLEN_T,MF_Pre,MF_Suf,TD,(File_Info.FILETYPE==FQ)? Head.Quality:Dummy_Q);
					TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr].p=UINT_MAX;
					if(TD.Compiled_Junctions_Ptr)
					{
						Real_Skip= -Skip_Bases;
						Skip_Bases=0;
						break;
					}
				}
			}
		}
		else if(Mismatch_Scan>=0)//revert to match..
		{
			Print_Matches(MF,MC,Head,SAM[0],READLEN_T,Genome_Offsets,Mismatch_Scan);
			continue;
		}

		if(TD.Compiled_Junctions_Ptr)
		{
			if(!Err)
			{
				int firstSignal = -2;
				int tempType = Classify_Hits(TD.Compiled_Junctions,firstSignal);
				int approvedPtr;
				approvedPtr = getBest(Head.Tag+Skip_Bases,READLEN_T,TD.Compiled_Junctions, selectedJunctions, true,(File_Info.FILETYPE==FQ)? Head.Quality+Skip_Bases:Dummy_Q);

				if(approvedPtr > 1) 
				{
					Hit_ID++;
					if(ONEMULTIHIT) approvedPtr=1;
					for(int i=0; i<approvedPtr; i++) 
					{
						Print_Hits(Head,TD.Compiled_Junctions,SAM[0],selectedJunctions[i],Hit_ID,Err,Genome_Offsets,true,READLEN_T,Real_Skip);
					}
				}	
				else if(approvedPtr) 
				{
					for(int i=0; i<approvedPtr; i++) 
					{
						Print_Hits(Head,TD.Compiled_Junctions,SAM[0],selectedJunctions[i],0,Err,Genome_Offsets,false,READLEN_T,Real_Skip);
					}
				}
				else
				{
					SAM[0] << SAM_Read.SAM_Line;
					//Print_Unmapped(Head,File_Info.STRINGLENGTH,SAM[0]);
				}
			} 
			else
			{
				SAM[0] << SAM_Read.SAM_Line;
				//Print_Unmapped(Head,File_Info.STRINGLENGTH,SAM[0]);
			}

		}
		else if(TD.Partial_Junctions_Ptr)
		{
			int Type=0;
			if(TD.Partial_Junctions_Ptr==1 && TD.Partial_Junctions[0].Junc_Count==1 && (Type=Canonical_Score(TD.Partial_Junctions[0].signal)))
			{
				assert(Type);
				selectedJunctions[0]=0;TD.Partial_Junctions[0].Type=Type;
				Print_Hits(Head,TD.Partial_Junctions,SAM[0],selectedJunctions[0],0,Err,Genome_Offsets,false,READLEN_T,Real_Skip);
			}
			else
			{
				SAM[0] << SAM_Read.SAM_Line;
				//Print_Unmapped(Head,File_Info.STRINGLENGTH,SAM[0]);
			}
		}
		else
		{
			SAM[0] << SAM_Read.SAM_Line;
			//Print_Unmapped(Head,File_Info.STRINGLENGTH,SAM[0]);
		}
	}
	Join_Tables(Genome_Offsets,Thread_ID);

}

void Scan_End_Junc(char* Current_Tag,int StringLength,Transcript_Data & TD,int & Err,char* Q)
{
	char Rev[StringLength],Rev_Bin[StringLength],RQ[MAXTAG];
	char Fwd[StringLength];
	for(int i=0;i<StringLength;i++)
	{
		Fwd[i]="ACGT"[Current_Tag[i]];
	}
	Convert_Reverse(Current_Tag,Rev,Rev_Bin,StringLength);
	Reverse_Quality(RQ,Q,StringLength);

	TD.Compiled_Junctions_Ptr=0;
	LEN L_Temp=TD.Generic_Hits.L;

	Do_End_Scan(Current_Tag,Fwd,StringLength,1,TD,Err,Q);//+ strand.. 
	Do_End_Scan(Rev_Bin,Rev,StringLength,0,TD,Err,RQ);//+ strand.. 
	
	TD.Generic_Hits.L=L_Temp;
	TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr].p=UINT_MAX;
}

const int MAX_PATTERN_MATCH=5;
void Do_End_Scan(char* Current_Tag,char* Current_Tag_ASCII,int StringLength,char Sign,Transcript_Data & TD,int & Err,char* Q) 
{
	if(Err) return;
	char* Middle_Tag=Current_Tag+MINX;
	char* Middle_Tag_ASCII=Current_Tag_ASCII+MINX;
	int Middle_StringLength=StringLength-2*MINX;
	char Org_Ten_MerR[TENMER+1];Org_Ten_MerR[TENMER]=0;
	char Org_Ten_MerL[TENMER+1];Org_Ten_MerL[TENMER]=0;
	
	for (int i=0;i<TENMER;i++)
	{
			Org_Ten_MerR[i]=Current_Tag_ASCII[i+StringLength-TENMER];
	}
	for (int i=0;i<TENMER;i++)
	{
			Org_Ten_MerL[i]=Current_Tag_ASCII[i];
	}

	Reset_MEMX(TD.Generic_Hits);TD.Generic_Hits.FMIndex=REVERSE;
	LEN L;L.IGNOREHEAD=0;
	Split_Read(Middle_StringLength,L);TD.Generic_Hits.L=L;
	TD.Generic_Hits.Current_Tag=Middle_Tag;
	
	int Mismatches=Scan(TD.Generic_Hits,2,L,fwfmi,revfmi,0,UINT_MAX);
	if(Mismatches<0) return;//No match

	char LeftX[100],Right[100];char* Left=LeftX+2;
	int Compiled_Junctions_Ptr=TD.Compiled_Junctions_Ptr;
	int Pattern_Match=0;
	
	for(int i=0;TD.Generic_Hits.Hit_Array[i].Start && !Err;i++)
	{
		SARange SA=TD.Generic_Hits.Hit_Array[i];
		assert(SA.Start && TD.Generic_Hits.Hits>0);
		if(SA.Start==SA.End)
		{
			if(!Get_Bases_ASCII(SA.Start-Middle_StringLength+RQFACTOR-MINX-2,MINX+2,LeftX) || !Get_Bases_ASCII(SA.Start+RQFACTOR,MINX+1,Right))
			{
				Err++;
				return;
			}
			Left[MINX+2]=0;//Right[MINX+1]=0;
			int Left_Mis=0,Right_Mis=0;
			float LScore=0,RScore=0;
			for(int j=0;j<MINX;j++)
			{
				if(Left[j]!=Current_Tag_ASCII[j]) 
				{
					LScore+=Penalty(Q[j]);
					Left_Mis++;
				}
				if(Right[j]!=Current_Tag_ASCII[j+Middle_StringLength+MINX]) 
				{
					Right_Mis++;
					RScore+=Penalty(Q[j]);
				}
			}

			//if(Left_Mis<=MIS_DENSITY)//Left anchors well..
			if(LScore<=QSUM_LIMIT)//Left anchors well..
			{
				char Right_Long[EXON_GEN_LEN+1];//Right_Long[EXON_GEN_LEN]=0;
				char* Ten_Mer=Right_Long;
				Ann_Info A1,A2;
				unsigned Loc1=SA.Start+RQFACTOR-Middle_StringLength;
				unsigned Loc2=SA.Start+RQFACTOR-Middle_StringLength+EXON_GEN_LEN+MINX;

				Location_To_Genome(Loc1,A1);
				Location_To_Genome(Loc2,A2);
				if (Loc1+StringLength > A1.Size||A1.ID!=A2.ID)//check for a Boundary Hit..
				{
					continue;
				}

				if(!Get_Bases_ASCII(SA.Start+RQFACTOR,EXON_GEN_LEN,Right_Long))
				{
					Err++;
					return;
				}
				while(Ten_Mer)
				{
					Ten_Mer=strstr(Ten_Mer,Org_Ten_MerR);
					if(Ten_Mer)//TODO:handle matches..
					{
						Pattern_Match++;
						if(Pattern_Match>MAX_PATTERN_MATCH ||(Ten_Mer-Right_Long)<80)
						{
							Err++;break;
						}
						Junction* junctions = extendX(Current_Tag_ASCII+MINX+Middle_StringLength-1,Right,Ten_Mer-(MINX-TENMER),SA.Start+RQFACTOR-1,SA.Start+RQFACTOR+(Ten_Mer-Right_Long),Sign,MINX+Middle_StringLength);
						for(int i=0;junctions[i].p!=UINT_MAX;i++)
						{
							junctions[i].Label=Compiled_Junctions_Ptr;junctions[i].Junc_Count=1;junctions[i].ID=INT_MAX-1;junctions[i].Sign=Sign;

							if(TD.Compiled_Junctions_Ptr<MAX_JUNCS_ALLOWED)
							{
								TD.Compiled_Junctions[Compiled_Junctions_Ptr++]=junctions[i];
							}
							else
							{
								Err++;break;
							}
						}
						delete [] junctions;
						Ten_Mer+=TENMER;
					}

				}

			}

			//if(Right_Mis<=MIS_DENSITY)//Right anchors well..
			if(RScore<=QSUM_LIMIT)//Right anchors well..
			{
				char Left_Long[EXON_GEN_LEN+1];//Right_Long[EXON_GEN_LEN]=0;
				char* Ten_Mer=Left_Long;
				Ann_Info A1,A2;
				unsigned Loc1=SA.Start+RQFACTOR-Middle_StringLength;
				unsigned Loc2=SA.Start+RQFACTOR-Middle_StringLength-EXON_GEN_LEN;

				Location_To_Genome(Loc1,A1);
				Location_To_Genome(Loc2,A2);
				if (Loc1+StringLength > A1.Size||A1.ID!=A2.ID)//check for a Boundary Hit..
				{
					continue;
				}
				if(!Get_Bases_ASCII(SA.Start+RQFACTOR-Middle_StringLength-EXON_GEN_LEN,EXON_GEN_LEN,Left_Long))
				{
					Err++;
					return;
				}
				while(Ten_Mer)
				{
					Ten_Mer=strstr(Ten_Mer,Org_Ten_MerL);
					if(Ten_Mer)//TODO:handle matches..
					{
						Pattern_Match++;
						if(Pattern_Match>MAX_PATTERN_MATCH ||(EXON_GEN_LEN-(Ten_Mer-Left_Long)<80))
						{
							Err++;break;
						}

						Junction* junctions = extendX(Current_Tag_ASCII+TENMER-1,Ten_Mer+TENMER,Left+TENMER,SA.Start+RQFACTOR-Middle_StringLength-EXON_GEN_LEN+(Ten_Mer-Left_Long)+TENMER-1,SA.Start+RQFACTOR-Middle_StringLength,Sign,TENMER);
						for(int i=0;junctions[i].p!=UINT_MAX;i++)
						{
							junctions[i].Label=Compiled_Junctions_Ptr;junctions[i].Junc_Count=1;junctions[i].ID=INT_MAX-1;junctions[i].Sign=Sign;
							if(TD.Compiled_Junctions_Ptr<MAX_JUNCS_ALLOWED)
							{
								TD.Compiled_Junctions[Compiled_Junctions_Ptr++]=junctions[i];
							}
							else
							{
								Err++;break;
							}
						}
						delete [] junctions;
						Ten_Mer+=TENMER;
					}

				}

			}
		}
	}

	TD.Compiled_Junctions_Ptr=Compiled_Junctions_Ptr;
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

int getBest(char* Current_Tag,int StringLength,Junction * Final_Juncs, int * approved, bool print,char* Q) {
	float max = -100000;
	int ptr = 0;
	int count = 0;
	char RC_Read[MAXDES],RC_Bin[MAXDES],RQ[MAXTAG];
	char Cat[MAXDES];
	int Mismatch_Pos[MAXDES];
	Convert_Reverse(Current_Tag,RC_Read,RC_Bin,StringLength);
	Reverse_Quality(RQ,Q,StringLength);

//Calculate true mismatches..
	int Junc_Count=Final_Juncs[0].Junc_Count;
	for(int i =0; Final_Juncs[i].p != UINT_MAX;i+=(Junc_Count ? Junc_Count:1),\
						   Junc_Count=Final_Juncs[i].Junc_Count) 
	{
		bool Full_String=false;
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
			Full_String=true;
			Get_Bases_ASCII(Final_Juncs[i].p, StringLength,Cat);
		}

		char* Read;char* Qual;
		Final_Juncs[i].Mismatches=0;
		if(Final_Juncs[i].Sign)
		{
			Read=Current_Tag;Qual=Q;
		}
		else
		{
			Read=RC_Bin;Qual=RQ;
		}

		int Pos=0;float PScore=0;
		for(int j=0;j<StringLength;j++)
		{
			if(Cat[j]!="ACGT"[Read[j]]) 
			{
				Final_Juncs[i].Mismatches++;
				PScore+=Penalty(Qual[j]);
				Mismatch_Pos[Pos++]=j;
			}
			Cat[j]=0;
		}

		if(Full_String && PScore<QSUM_LIMIT)
		{
			return 0;
		}
		if(PScore>QSUM_LIMIT)
		{
			Final_Juncs[i].Mismatches==100;
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

		if(!Final_Juncs[i].Mismatches==100) continue; 
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
				//Junc_Score+=Final_Juncs[i+j].Type-Junc_Count*2;
				if(dist<(Final_Juncs[i+j].q-Final_Juncs[i+j].p)) //Pick Min dist
				{
					dist=Final_Juncs[i+j].q-Final_Juncs[i+j].p;
				}
				//dist+=Final_Juncs[i+j].q-Final_Juncs[i+j].p;
			}
			Junc_Score-=Junc_Count*2;
			//dist=dist/Junc_Count;

			assert(dist>0);
			if(dist>30000) Junc_Score+=1;
			else if(dist>5000) Junc_Score+=2;
			else if(dist>2000) Junc_Score+=3;
			else if(dist>=100) Junc_Score+=4;
			else if(dist>60) Junc_Score+=1;
		}
		int tempScore = -3*Final_Juncs[i].Mismatches+Junc_Score;//+Final_Juncs[i].score;
		if(Final_Juncs[i].Junc_Count<=3 && tempScore >= max)
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

	if (!MIS_DENSITY || Err)
	{
		Final_Juncs[Final_Juncs_Ptr].p=UINT_MAX;
		return Final_Juncs_Ptr;
	}
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

bool Align(char* Source,int StringLength,unsigned Dest_Loc,SARange & R,int Actual_Length,int Read_Skip,Transcript_Data & TD,char* Q)
{


	if(Dest_Loc!=UINT_MAX)//extend in genome..
	{
		char Dest[StringLength];
		Get_Bases(Dest_Loc,StringLength,Dest);
		int Mis_Count=0;
		float PScore=0;
		//for(int i=0;i<StringLength && Mis_Count<=1/*Max mis allowed*/;i++)
		for(int i=0;i<StringLength && PScore<=QSUM_LIMIT;i++)
		{
			assert(Dest[i]<=4 && Dest[i]>=0 && Source[i]<=4 && Source[i]>=0);
			if(Dest[i]!=Source[i])
			{
				Mis_Count++;
				PScore+=Penalty(Q[i]);
			}
		}
		//if (Mis_Count>MIS_DENSITY+1)
		if(PScore>QSUM_LIMIT)
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
		if(DEEPSCAN && !Read_Skip)
		{
			TD.Generic_Hits.Current_Tag=Source;TD.Generic_Hits.Hit_Array_Ptr=0;TD.Generic_Hits.Hits=0;
			int Last_Mis=Deep_Scan(TD.Generic_Hits,MIS_DENSITY+1,TD.Generic_Hits.L,fwfmi,revfmi,UINT_MAX);
			if(TD.Generic_Hits.Hit_Array_Ptr) return true;
			else return false;
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
}

inline void B2C(char *S,char *D,int StringLength)
{
	*(D+StringLength)=0;
	for(int i=0;i<StringLength;i++,S++,D++)
		*D="ACGT"[*S];
}

int Seek_Junc(char* S,SARange R,int Read_Skip,int Junc_Count,int Mis_In_Junc_Count,unsigned Last_Exon,int StringLength,int Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int & Inspected_Pairs,int & Err,bool Sign,Transcript_Data & TD,char* Q)
{
	assert(Trans_Array_Ptr<MAX_JUNCS_IN_TRANSCRIPT && TD.Compiled_Junctions_Ptr>=0);assert(Read_Skip>=0 && StringLength >=0 && StringLength<=READLEN && Inspected_Pairs>=0);
	assert(Read_Skip>=MINX || Read_Skip==0);
	if(Read_Skip>StringLength-MINX)//End of read..
	{
		if((Read_Skip==StringLength) || Align(S+Read_Skip,StringLength-Read_Skip,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip,R,StringLength,Read_Skip,TD,Q+Read_Skip))//End of the read can be matched..
		{
			if(!Junc_Count)//Match..
			{
				TD.Max_Junc_Found=Junc_Count;
				Junction J;J.Sign=(Sign)? 1:0;J.Mismatches=R.Mismatches;J.Junc_Count=0;
				if (Last_Exon!=UINT_MAX)//one hit.. 
				{
					J.p=Last_Exon+1;
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

	if(Align(S+Read_Skip,MINX,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip,R,Read_Skip+MINX,Read_Skip,TD,Q+Read_Skip))//Can exact extension be done..
	{
		//if(TD.Generic_Hits.Hit_Array_Ptr>100) {TD.Generic_Hits.Hit_Array_Ptr=100;}//Err++;return 0;}
		std::vector <SARange> Hits(TD.Generic_Hits.Hit_Array_Ptr);
		int Ptr=0;
		for (int i=0;DEEPSCAN && Last_Exon==UINT_MAX && i<TD.Generic_Hits.Hit_Array_Ptr;i++)//Save hits in a vector..
		{
			assert(TD.Generic_Hits.Hit_Array[i].Start);
			SARange SA=TD.Generic_Hits.Hit_Array[i];
			float Score=0;
			for(int j=0;j<SA.Mismatches;j++)
			{
				Score+=Penalty(Q[j+Read_Skip]);
			}
			if(Score<32)
			{
				Hits[Ptr++]=SA;
			}

		}
		if(!Ptr)
		{
			for (int i=0;i<TD.Generic_Hits.Hit_Array_Ptr;i++)//Save hits in a vector..
			{
				assert(TD.Generic_Hits.Hit_Array[i].Start);
				Hits[i]=TD.Generic_Hits.Hit_Array[i];
			}
		}
		else
			TD.Generic_Hits.Hit_Array_Ptr=Ptr;

		unsigned Last_ExonT;int Hit_Count=TD.Generic_Hits.Hit_Array_Ptr-1;
		for (int i=0;i<=Hit_Count;i++)//Save hits in a vector..
		{
			SARange SA=Hits[i];
			Last_ExonT=Last_Exon;
			if(SA.Start==SA.End && Last_Exon==UINT_MAX) 
			{
				Last_ExonT=SA.Start-Read_Skip;
			}
			Seek_Junc(S,SA,Read_Skip+MINX,Junc_Count,Mis_In_Junc_Count,Last_ExonT,StringLength,Trans_Array_Ptr,Pre,Suf,Inspected_Pairs,Err,Sign,TD,Q);
			if(Err) return 0;
			if (Inspected_Pairs >= MAX_INSPECTED_PAIRS) {Err++;return 0;}
		}
	}

//Junction check routines
	if(Read_Skip>=MINX)
	{
		if(Exons_Too_Small(Trans_Array_Ptr,TD)) return 0;
		//Case of a single junction between S[a,a+2l-1]
		if (Junc_Count<=TD.Max_Junc_Found) //This will add one junction..
		{
			//the junction is between S[Read_Skip,Read_Skip+18]
			//so anchor S[Read_Skip-18,Read_Skip]
			//Penguin S[Read_Skip-18,Read_Skip] ---- S[Read_Skip+18,Read_Skip+36] and extend..
			Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX,TD,Q+Read_Skip-MINX);
			//Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX-6,TD,Q+Read_Skip-MINX);
			//Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX-12,TD,Q+Read_Skip-MINX);
			//Enum_Single_Junctions(S,S+Read_Skip-MINX,Read_Skip,StringLength,(Last_Exon==UINT_MAX) ? UINT_MAX:Last_Exon+Read_Skip-MINX,Inspected_Pairs,R,Junc_Count,Mis_In_Junc_Count,Trans_Array_Ptr,Pre,Suf,Read_Skip-MINX,Err,Sign,3*MINX-15,TD,Q+Read_Skip-MINX);
		}
		if(!TD.Compiled_Junctions_Ptr)
		{
			if(TD.Ext_Array_Ptr<MAX_JUNCS_ALLOWED)
			{
				TD.Ext_Array[TD.Ext_Array_Ptr].Read_Skip=Read_Skip;
				TD.Ext_Array[TD.Ext_Array_Ptr].Sign=Sign;
				TD.Ext_Array[TD.Ext_Array_Ptr++].SA=R;
			}
		}
	}

	return Trans_Array_Ptr;

}

void Enum_Single_Junctions(char* Org_Read,char* Converted_Read,int Read_Skip,int StringLength, unsigned Anchor,int & Inspected_Pairs,SARange & R,int Junc_Count,int Mis_In_Junc_Count,int & Trans_Array_Ptr,MEMX & Pre,MEMX & Suf,int Skip,int & Err,int Sign,int Fragment,Transcript_Data & TD,char* Q)
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
					bool Existing_Junc=true;
					int Label=TD.Compiled_Junctions[j].Label;
					if(j-Trans_Array_Ptr<0) Existing_Junc=false;
					for(int k=1;k<=Trans_Array_Ptr && Existing_Junc;k++)//Scan if prev. junc the same..
					{
						if(Label!=TD.Compiled_Junctions[j-k].Label)
							Existing_Junc=false;
						else if(TD.Compiled_Junctions[j-k].p!=TD.Trans_Array[Trans_Array_Ptr-k].p || TD.Compiled_Junctions[j-k].q!=TD.Trans_Array[Trans_Array_Ptr-k].q)
						{
							Existing_Junc=false;
						}
					}
					if(Existing_Junc)
					{
						if(j!=Trans_Array_Ptr)
						{
							if(Label==TD.Compiled_Junctions[j-Trans_Array_Ptr-1].Label)//Actually there is another extra junc.
							{
								continue;
							}
						}
						Junc_Already_Found=true;break;
					}
				}
			}
			if(Junc_Already_Found) continue;
			if(!PRINT_NON_CANON && !Canonical_Score(Final_Juncs[i].signal)) 
				continue;
			TD.Trans_Array[Trans_Array_Ptr]=Final_Juncs[i];//.x=p;Trans_Array[Trans_Array_Ptr].y=q;
			Seek_Junc(Converted_Read+r-Skip,R,0,Junc_Count+1,Mis_In_Junc_Count,q+1,StringLength-r,Trans_Array_Ptr+1,Pre,Suf,Inspected_Pairs,Err,Sign,TD,Q+r-Skip);
			if(Err) break;
			if (Inspected_Pairs >= MAX_INSPECTED_PAIRS) {Err++;break;}
		}
	}
	else Err++;
	delete [] Pairs;
}

int Seek_Single_Strand(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,int Sign,Transcript_Data & TD,char* Q)
{
	SARange SA;
	SA.Start=0;SA.End=revfmi->textLength;//SA.FMIndex=REVERSE;
	SA.Level=1; SA.Skip=0;SA.Mismatches=0;SA.Mismatch_Char=0;Reset_MEMX(TD.Generic_Hits);TD.Generic_Hits.Hit_Array_Ptr=0;TD.Generic_Hits.Current_Tag=Current_Tag;

	int Inspected_Pairs=0;
	int Err=0;
	Seek_Junc(Current_Tag,SA,0,0,0,UINT_MAX,StringLength,0,MF_Pre,MF_Suf,Inspected_Pairs,Err,Sign,TD,Q);
	return Err;
}

int Seek_All_Junc(char* Current_Tag,int StringLength,MEMX & MF_Pre,MEMX & MF_Suf,Transcript_Data & TD,char* Q)
{
	int Err= Seek_Single_Strand(Current_Tag,StringLength,MF_Pre,MF_Suf,1/*Plus*/,TD,Q);
	if (Err) return Err;
	char RC_Read[MAXDES],RC_Bin[MAXDES],RQ[MAXDES];
	Convert_Reverse(Current_Tag,RC_Read,RC_Bin,StringLength);
	Reverse_Quality(RQ,Q,StringLength);
	Err+= Seek_Single_Strand(RC_Bin,StringLength,MF_Pre,MF_Suf,0/*Minus*/,TD,RQ);
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

float Pr(float Q)
{
	assert((int)Q>=0 && (int)Q<POWLIMIT-1);
	return(1-POW10[(int)Q]);
	//printf("table: %f\tlib: %f\n",POW10[(int)Q],1-pow(10,-Q/10));
	//return (1-pow(10,-Q/10));
}

float Pow10(float Q)
{
	static float Max=0;
	//if(Max<Q) {Max=Q;printf("Newmax: %f\n",Max);}
	assert((int)Q>=0);
	if((int)Q<POWLIMIT-1)
		return POW10[(int)Q];
	else
		return pow(10,-float(Q)/10);
}

void Build_Pow10()
{
	for(int Q=0;Q<POWLIMIT;Q++)
	{
		POW10[Q]=(pow(10,-float(Q)/10));
	}
}


float Penalty(char Q)
{
	float Q_Value=Q-QUALITYCONVERSIONFACTOR;//make quality to integer..
	//Q_Value=std::min(30.0f,Q_Value);
	assert(Q_Value>=0);// && Q_Value<=40);
	float Penalty= -10*log10((1-Pr(Q_Value))/3);
	Penalty=std::min(QLIMIT_FLOAT,Penalty);
	return Penalty;
}

void Reverse_Quality(char* Dest,char* Q,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=Q[i];Dest++;
	}
	*Dest=0;
}

bool Mismatch_Hit_Nice(int Mismatch_Scan,MEMX & MF,MEMX & MC,char* Q,int StringLength)
{
	char RQ[MAXTAG];
	if(Mismatch_Scan== -1)//No mismatch hits..
	{
		return false;
	}
	SARange SA;float PScore=0;
	if(MF.Hits+MC.Hits>1)//Many mismatch hits..
		return true;
	if(MF.Hits)
	{
		SA=MF.Hit_Array[0];
		for (int i=0;i<SA.Mismatches;i++)
		{
			PScore+=Penalty(Q[SA.Mismatch_Pos[i]]);
		}
	}
	else
	{
		Reverse_Quality(RQ,Q,StringLength);
		SA=MC.Hit_Array[0];
		for (int i=0;i<SA.Mismatches;i++)
		{
			PScore+=Penalty(Q[SA.Mismatch_Pos[i]]);
		}
	}
	if(PScore<30) 
		return true;
	else
		return false;
}

/*void Scan_Partial_Read(char* Current_Tag,int StringLength,MEMX & Pre, MEMX & Suf,Transcript_Data & TD,char* Q)
{
	char Rev[StringLength],Rev_Bin[StringLength],RQ[MAXTAG];
	char Fwd[StringLength];
	for(int i=0;i<StringLength;i++)
	{
		Fwd[i]="ACGT"[Current_Tag[i]];
	}
	Convert_Reverse(Current_Tag,Rev,Rev_Bin,StringLength);
	Reverse_Quality(RQ,Q,StringLength);

	TD.Compiled_Junctions_Ptr=0;
	int Err=0;
	Partial_Scan(Current_Tag+StringLength-50,Fwd+StringLength-50,50,Pre,Suf,1,TD,Err,Q);//+ strand.. 
	Partial_Scan(Rev_Bin+StringLength-50,Rev+StringLength-50,50,Pre,Suf,0,TD,Err,RQ);//+ strand.. 
	TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr].p=UINT_MAX;
}

void Partial_Scan(char* Current_Tag,char* Current_Tag_ASCII,int StringLength,MEMX & Pre,MEMX & Suf,char Sign,Transcript_Data & TD,int & Err,char* Q) 
{
	if(Err) return;
	char* Head_Tag=Current_Tag;
	int Head_StringLength=StringLength/2;
	char Org_18_MerR[MINX+1];Org_18_MerR[MINX]=0;
	
	for (int i=0;i<MINX;i++)
	{
			Org_18_MerR[i]=Current_Tag_ASCII[i+StringLength-MINX];
	}

	Reset_MEMX(TD.Generic_Hits);TD.Generic_Hits.FMIndex=REVERSE;
	LEN L,L_Temp;L_Temp=TD.Generic_Hits.L;L.IGNOREHEAD=0;
	Split_Read(Head_StringLength,L);TD.Generic_Hits.L=L;
	TD.Generic_Hits.Current_Tag=Head_Tag;

//scan first half..	
	int Mismatches=Scan(TD.Generic_Hits,2,L,fwfmi,revfmi,0,UINT_MAX);
	if(Mismatches<0) return;

	char LeftX[100],Right[100];char* Left=LeftX+2;
	int Compiled_Junctions_Ptr=TD.Compiled_Junctions_Ptr;
	int Pattern_Match=0,Label=0;
	Err=0;
	LEN L18;
	L18.IGNOREHEAD=0;Split_Read(MINX,L18);//we are scanning MINX-mers...
	PAIR* Pairs= new PAIR[MAX_HITS_TO_STORE+10];if(!Pairs) {cout << "Seek_Junc(): Error allocatind memory..\n";exit(100);}
	Junction Final_Juncs[MAX_JUNCS_TO_STORE+1];Final_Juncs[0].p=UINT_MAX;
	int Inspected_Pairs=0;
	
	for(int i=0;TD.Generic_Hits.Hit_Array[i].Start && !Err;i++)
	{
		SARange SA=TD.Generic_Hits.Hit_Array[i];
		assert(SA.Start && TD.Generic_Hits.Hits>0);
		if(SA.Start==SA.End)
		{
			char Read[200];
			char Read_Head[200];
			B2C(Current_Tag,Read,StringLength);B2C(Current_Tag,Read_Head,MINX);
			Inspected_Pairs+=Find_Single_Junc(Read,Read_Head,Current_Tag,Pre,Suf,StringLength,L18,Pairs,Final_Juncs,Err,Label,'+',SA.Start,SA,1);
		}
	}

	delete [] Pairs;
	TD.Generic_Hits.L=L_Temp;
	TD.Compiled_Junctions_Ptr=Compiled_Junctions_Ptr;
}*/

const int MAX_LOC=5;
const int SEED_LEN=14;
void Scan_Right_Read(char* Current_Tag,int StringLength,MEMX & Pre, MEMX & Suf,Transcript_Data & TD,char* Q)
{
	char Rev[StringLength],Rev_Bin[StringLength],RQ[MAXTAG];
	char Fwd[StringLength];
	for(int i=0;i<StringLength;i++)
	{
		Fwd[i]="ACGT"[Current_Tag[i]];
	}
	Convert_Reverse(Current_Tag,Rev,Rev_Bin,StringLength);
	Reverse_Quality(RQ,Q,StringLength);

	TD.Compiled_Junctions_Ptr=0;
	int Err=0;
	Extension_Scan(Current_Tag,Fwd,StringLength,Pre,Suf,1,TD,Err,Q);//+ strand.. 
	if(!Err)
		Extension_Scan(Rev_Bin,Rev,StringLength,Pre,Suf,0,TD,Err,RQ);//+ strand.. 
	TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr].p=UINT_MAX;
}

void Extension_Scan(char* Current_Tag,char* Current_Tag_ASCII,int StringLength,MEMX & Pre,MEMX & Suf,char Sign,Transcript_Data & TD,int & Err,char* Q) 
{
	unsigned Loc[MAX_LOC+1];
	char GTAG[StringLength],ATAC[StringLength],GCAG[StringLength];
	int GTAG_Ptr=0,GCAG_Ptr=0,ATAC_Ptr=0;
	for(int r=0;r<TD.Ext_Array_Ptr;r++)
	{
		if(TD.Ext_Array[r].Sign !=Sign)
			continue;
		SARange SA=TD.Ext_Array[r].SA;
		assert(SA.FMIndex==REVERSE);
		int Offset=TD.Ext_Array[r].Read_Skip;
		int Loc_Ptr=0;
		if(SA.Start==SA.End)
		{
			Loc[Loc_Ptr++]=SA.Start;
		}
		else
		{
			for(Loc_Ptr=0;SA.Start+Loc_Ptr<=SA.End;Loc_Ptr++)
			{
				Loc[Loc_Ptr]=revfmi->textLength-Offset-BWTSaValue(revfmi,SA.Start+Loc_Ptr);
				if(Loc_Ptr==MAX_LOC)
				{
					Err++;
					return;
				}
			}
		}

		char RightX[EXON_GEN_LEN+1];
		int GT[StringLength],CT[StringLength],AT[StringLength],GC[StringLength];
		for(int i=0;i<Loc_Ptr;i++)
		{
			int GT_Ptr=0,CT_Ptr=0,AT_Ptr=0,GC_Ptr=0;
			if(!Get_Bases_ASCII(Loc[i],StringLength-SEED_LEN,RightX))
			{
				Err++;return;
			}

			RightX[StringLength-SEED_LEN]=0;
			char* Scan_Read=RightX+Offset-1;
			char* Last_Offset;
			float Match_Penalty=0;
			for(int Last_Offset=0;RightX+Last_Offset==Scan_Read;Last_Offset++)
			{
				if(*(RightX+Last_Offset)==Current_Tag_ASCII[Last_Offset])
				{
					Match_Penalty+=Penalty(Q[Last_Offset]);
				}
			}

			while((Scan_Read=strpbrk(Scan_Read+1,"TC")) && (Match_Penalty<QSUM_LIMIT))
			{
				for(int Last_Offset=0;RightX+Last_Offset==Scan_Read;Last_Offset++)
				{
					if(*(RightX+Last_Offset)==Current_Tag_ASCII[Last_Offset])
					{
						Match_Penalty+=Penalty(Q[Last_Offset]);
					}
				}
				if(Match_Penalty>=QSUM_LIMIT)
					continue;

				char Prec=*(Scan_Read-1);
				assert(Prec=='A'||Prec=='C'||Prec=='G'||Prec=='T');
				if(*Scan_Read=='T')
				{
					switch(Prec)
					{
						case 'A':
							AT[AT_Ptr++]=Scan_Read-RightX-1;
							break;
						case 'C':
							CT[CT_Ptr++]=Scan_Read-RightX-1;
							break;
						case 'G':
							GT[GT_Ptr++]=Scan_Read-RightX-1;
					}
				}
				else
				{
					if(Prec=='G')
					{
						GC[GC_Ptr]=Scan_Read-RightX-1;
					}
				}
			}

			if(GC_Ptr+GT_Ptr+CT_Ptr+AT_Ptr)//Motifs present..
			{
				if(!Get_Bases_ASCII(Loc[i]+StringLength-SEED_LEN,EXON_GEN_LEN-StringLength+SEED_LEN,RightX+StringLength-SEED_LEN))//load rest..
				{
					Err++;return;
				}
				RightX[EXON_GEN_LEN]=0;
				char Signal[5];Signal[4]=0;

				if(GT_Ptr)
				{
					Signal[0]='G';Signal[1]='T';
					Seek_Residue(GT_Ptr,GT,Current_Tag_ASCII,StringLength,RightX,Q,Signal,Loc[i],TD,Sign);
				}
				if(CT_Ptr)
				{
					Signal[0]='C';Signal[1]='T';
					Seek_Residue(CT_Ptr,CT,Current_Tag_ASCII,StringLength,RightX,Q,Signal,Loc[i],TD,Sign);
				}


			}

			/*Ann_Info A1;
			Location_To_Genome(Loc[r],A1);
			cout<<A1.Name<<"\t"<<Loc[r]<<endl;*/
			
		}
		
	}
}

void Seek_Residue(int & GT_Ptr,int *GT,char* Current_Tag_ASCII,int StringLength,char* RightX,char* Q,char *Signal,unsigned Loc,Transcript_Data & TD,char Sign)
{
	char Seed[SEED_LEN+1];Seed[SEED_LEN]=0;
	char* Locate;
	for(int i=0;i<GT_Ptr;i++)
	{
		memcpy(Seed,Current_Tag_ASCII+GT[i],SEED_LEN);
		Locate=RightX+StringLength;
		while(Locate=strstr(Locate+1,Seed))
		{
			float Score=0;
			int Mis=0,j=0;
			for(int k=GT[i];k<StringLength && Score<QSUM_LIMIT;k++,j++)
			{
				if(Locate[j]!=Current_Tag_ASCII[k])
				{
					Mis++;
					Score+=Penalty(Q[k]);
				}
			}
			if(Score<QSUM_LIMIT)
			{
				Signal[2]=*(Locate-2);
				Signal[3]=*(Locate-1); 
				if(int Type=Canonical_Score(Signal))
				{
					Junction J;
					J.Label=TD.Compiled_Junctions_Ptr;
					J.Junc_Count=1;J.ID=INT_MAX-1;J.Sign=Sign;J.Type=Type;
					J.p = Loc + StringLength-j;
					J.q = Loc + Locate-RightX-1;
					J.r = StringLength-j;
					strcpy(J.signal,Signal);
					if(TD.Compiled_Junctions_Ptr<MAX_JUNCS_ALLOWED)
					{
						TD.Compiled_Junctions[TD.Compiled_Junctions_Ptr++]=J;
					}
				}
			}
		}
	}
}


void Print_Matches(MEMX & MF,MEMX & MC,READ & Head,ofstream & SAM,int StringLength,Offset_Record *Genome_Offsets,int Mismatches)
{
	unsigned Loc,Conversion_Factor=revfmi->textLength-StringLength+1;
	if(MF.Hit_Array_Ptr) MF.Hit_Array_Ptr--;if(MC.Hit_Array_Ptr) MC.Hit_Array_Ptr--;
	bool Multi_Hits=(1==MF.Hit_Array_Ptr+MC.Hit_Array_Ptr)? false:true;
	Junction J;
	J.Mismatches=Mismatches;
	for(int i=0;i<MF.Hit_Array_Ptr;i++)
	{
		SARange SA=MF.Hit_Array[i];assert(SA.Start);
		if(SA.Start!=SA.End)
		{
			for(unsigned j=SA.Start;j<=SA.End;j++)
			{
				Loc=Conversion_Factor-BWTSaValue(revfmi,j);
				assert(Loc);
				J.p=Loc;J.r=0;J.q=0;J.Sign=1;
				Print_Hits(Head,&J,SAM,0,0,0,Genome_Offsets,Multi_Hits,READLEN,0);
				if(ONEMULTIHIT)
					return;
			}
		}
		else
		{
			Loc=SA.Start+RQFACTOR-StringLength+1;
			assert(Loc);
			J.p=Loc;J.q=J.r=0;J.Sign=1;
			Print_Hits(Head,&J,SAM,0,0,0,Genome_Offsets,Multi_Hits,READLEN,0);
			if(ONEMULTIHIT)
				return;
		}
	}	
	for(int i=0;i<MC.Hit_Array_Ptr;i++)
	{
		SARange SA=MC.Hit_Array[i];assert(SA.Start);
		if(SA.Start!=SA.End)
		{
			for(unsigned j=SA.Start;j<=SA.End;j++)
			{
				Loc=Conversion_Factor-BWTSaValue(revfmi,j);
				assert(Loc);
				J.p=Loc;J.r=0;J.q=0;J.Sign=0;
				Print_Hits(Head,&J,SAM,0,0,0,Genome_Offsets,Multi_Hits,READLEN,0);
				if(ONEMULTIHIT)
					return;
			}
		}
		else
		{
			Loc=SA.Start+RQFACTOR-StringLength+1;
			assert(Loc);
			J.p=Loc;J.q=J.r=0;J.Sign=0;
			Print_Hits(Head,&J,SAM,0,0,0,Genome_Offsets,Multi_Hits,READLEN,0);
			if(ONEMULTIHIT)
				return;
		}
	}	
}

void Print_Unmapped(READ & Head,int StringLength,ofstream & MISHIT)
{
	Nullify_String(Head.Description);
	Head.Quality[StringLength]=0;
	Head.Tag_Copy[StringLength]=0;
	MISHIT <<Head.Description+1<<"\t4\t*\t0\t0\t*\t*\t0\t0\t"<<Head.Tag_Copy<<"\t"<<Head.Quality<<endl;
}

bool Messy_CIGAR(char *Cig,int StringLength)
{
	char* Int_Start=Cig;
	int Len=0,Max_Len=0;
	for(;*Cig;Cig++)
	{
		if(*Cig=='M')
		{
			*Cig=0;
			Len=atoi(Int_Start);*Cig='M';
			if(Max_Len<Len)
			{
				Max_Len=Len;
			}
			Int_Start=Cig+1;
		}
		else if(!isdigit(*Cig))
		{
			Int_Start=Cig+1;
		}
	}
	assert(StringLength>=Max_Len);
	if(!Max_Len)
		return true;
	if(StringLength-Max_Len>=5)
		return true;
	else
		return false;
}
