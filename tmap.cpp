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
//#include "Hash.h"
#include "rqindex.h"
#include <fstream>
#include <iostream>
//#include "Print.h"
//#include "init.h"
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
#include "qcalc.h"
//#include "Cmdline.h"
#include "Indexes.h"
#include "file.h"
//#include "extend.h"
extern const int POWLIMIT=300;
extern const int QUALITYCONVERSIONFACTOR=33;

float QLIMIT_FLOAT=30.0f;
float QLIMIT= QLIMIT_FLOAT; 
char* LOG_SUCCESS_FILE=NULL;
FILE *Log_SFile;

Parameters CL;	
int THREAD;
float POW10[POWLIMIT];
Index_Info Genome_Files;
int MAX_MISMATCHES= -1;
MMPool *mmPool;
BWT *fwfmi,*revfmi;
unsigned CONVERSION_FACTOR;
unsigned Conversion_Factor,SOURCELENGTH;
map <unsigned, Ann_Info> Annotations;
FILETYPE File_Info;

int READLEN,ORG_STRINGLENGTH;
char Char_To_CodeC[256];
char Char_To_Code[256];

gzFile Input_File,Mate_File;
//---------------RQ stuff -------------------------
unsigned RQ_Hits;
SA* SA_Index;
int* SA_Blocks;
char COMPRESS;
int INDEX_RESOLUTION=30000;
int EXONGAP;
int RESIDUE=5;

void Parse_Command_line(int argc, char* argv[],Index_Info & Ind,Parameters & CL);
void Load_FM_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool);
void InitX(BWT *revfmi,unsigned & SOURCELENGTH,gzFile & Input_File,gzFile & Mate_File,FILETYPE & File_Info,Parameters & CL,Index_Info & Genome_Files,int & MAX_MISMATCHES);
void Convert_Reverse(char* Read,char * RC_Read,char* RC_Bin,int StringLength);
bool  Progress_Bar(Parameters & CL,unsigned & Number_of_Tags,unsigned & Progress,unsigned & Tag_Count,FILETYPE & File_Info);
void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T);
int Map_Read(MEMX & MF,MEMX & MC,int MAX_MISMATCHES, LEN & L,BWT* fwfmi,BWT* revfmi,int Next_Mis,int Max_Hits);
void *Map(void *T);
void Set_Affinity();
void Process_Hits(MEMX & MF,MEMX & MC,int StringLength,ofstream & Out_File,ofstream & Mishit_File,READ & Head);
void Open_Outputs(ofstream & SAM,string filename);
bool Call_Junc(char* Junc,unsigned Pos,int StringLength,ofstream & Out_File,string Des);

int main(int argc, char* argv[])
{

	time_t Start_Time,End_Time,Maptime;
	time(&Start_Time);

	Build_Pow10();
	Parse_Command_line(argc,argv,Genome_Files,CL);
	Load_FM_Indexes(Genome_Files,fwfmi,revfmi,mmPool);
	InitX(revfmi,SOURCELENGTH,Input_File,Mate_File,File_Info,CL,Genome_Files,MAX_MISMATCHES);

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

	//UnLoad_Indexes(fwfmi,revfmi,mmPool,Range_Index);
	fprintf(stderr,"\r[++++++++100%%+++++++++]\n");//progress bar....
	time(&End_Time);
	Maptime=difftime(End_Time,Start_Time);

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
	unsigned Total_Hits=0,Tags_Processed=0,Tag_Count=0;
	unsigned Number_of_Tags=1000;
	unsigned Progress=0;
	unsigned Hit_ID=0;

	MEMLOOK MLook;MLook.Lookupsize=6;
	Build_Tables(fwfmi,revfmi,MLook);
	LEN L;L.IGNOREHEAD=0;
	Split_Read(File_Info.STRINGLENGTH,L);//we are scanning whole read...
//--------------------- Setup Data Structure for Batman Prefix ----------------------------------------
	MEMX MF,MC;
	Init_Batman(MF,L,MLook,MAX_MISMATCHES);
	Init_Batman(MC,L,MLook,MAX_MISMATCHES);
//--------------------- Setup Data Structure for Batman End----------------------------------------
	
	ofstream Mishit_File,Out_File;
	string Mishit_File_Name = "rejected",Out_File_Name="trans";
	Open_Outputs(Mishit_File, Mishit_File_Name+Str_Thread_ID+".fq");
	Open_Outputs(Out_File, Out_File_Name+Str_Thread_ID+".sam");


	fprintf(stderr,"======================]\r[");//progress bar....
	int Actual_Tag=0;
	READ Head,Tail;
	while (Read_Tag(Head,Tail,Input_File,Mate_File,File_Info))
	{
		if(Head.NCount>2) continue;
		if(Thread_ID==1 && !Progress_Bar(CL,Number_of_Tags,Progress,Tag_Count,File_Info)) break;
		if(CL.MAX_TAGS_TO_PROCESS && CL.MAX_TAGS_TO_PROCESS<Actual_Tag) break;

		int Label=0;
		Actual_Tag++;

		char Rev_Bin[MAXTAG],Rev[MAXTAG];
		Convert_Reverse(Head.Tag,Rev,Rev_Bin,File_Info.STRINGLENGTH);
		MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Current_Tag=Head.Tag;//setup read details to alignmentstructure..
		MC.Hits=0;MC.Hit_Array_Ptr=0;MC.Current_Tag=Rev_Bin;MC.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
		int Mismatch_Scan=Map_Read(MF,MC,MAX_MISMATCHES,L,fwfmi,revfmi,0,2);
		if(Mismatch_Scan>=0)
		{
			Process_Hits(MF,MC,File_Info.STRINGLENGTH,Out_File,Mishit_File,Head);
		}
	}
}

void Process_Hits(MEMX & MF,MEMX & MC,int StringLength,ofstream & Out_File,ofstream & Mishit_File,READ & Head)
{
	int Plus_Hits=0,Minu_Hits=0;
	if(MF.Hit_Array_Ptr+MC.Hit_Array_Ptr!=1) //Multiple hits..
		return;
	
	Ann_Info A;
	unsigned Loc;
	SARange SA;
	if(MF.Hit_Array_Ptr)
	{
		SA=MF.Hit_Array[0];
	}
	else
		SA=MC.Hit_Array[0];

	if(SA.Start==SA.End)
	{
		Loc=SA.Start;

	}
	else//Multiple Hits..
	{
		return;
		Loc=Conversion_Factor-BWTSaValue(revfmi,SA.Start);
	}

	Location_To_Genome(Loc,A);
	if (Loc+StringLength > A.Size)//check for a Boundary Hit..
	{
		//string Des=A.Name;
		//Out_File<<A.Name<<":"<<Loc<<endl;
		*strchr(Head.Description,'\n')=0;
		string Des=Head.Description;
		char Name[300];strcpy(Name,A.Name);
		Call_Junc(Name,Loc,StringLength,Out_File,Des);
	}

}

option Long_Options[]=
{
{"help",0,NULL,'h'},
{"query",1,NULL,'q'},
{0,0,0,0}
};

void Parse_Command_line(int argc, char* argv[],Index_Info & Ind,Parameters & CL)
{
	int Current_Option=0;
	char Short_Options[] ="jM:bhq:t:g:G:n:N:o:w:mprO::T:";//allowed options....
	char* This_Program = argv[0];//Current program name....
	char Help_String[]=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --query | -m \t\t map only\n"
" --query | -M \t\t Store <number> junctions..\n"
" --query | -p \t\t process only\n"
" --query | -q <filename>\t\t Query file(File of Tags)\n"
" --query | -r \t\t mark exons in refGene..\n"
" --query | -o <filename>\t\t junction file...\n"
" --query | -w <filename>\t\t Wiggle file...\n"
" --query | -O <filename>\t\t log read mapping info...\n"
" --query | -G <number>\t\t maximum gap between two introns...\n"
" --query | -n <number>\t\t number of mismatches in Initial Scan...\n"
" --query | -N <number>\t\t number of mismatches in splices...\n"
" --query | -b \t\t build files from refgene\n"
" --query | -j \t\t output all junctions...\n"
;

	if(argc == 1) {printf("%s \n",Help_String);exit(0);}
	char *Source=(char*)malloc(sizeof(char)*6500);//create space for file names...
	char *options, *value; 
	char* Name;int Last_Dash;char* Genome_Name;
	char MARKEX=FALSE;

	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 'T':
				THREAD=atoi(optarg);
				break;
			case 't':
				CL.MAX_TAGS_TO_PROCESS=atoi(optarg);
				break;
			case 'q':
				if(!CL.Patternfile_Count){CL.PATTERNFILE=optarg;}
				else CL.PATTERNFILE1=optarg;
				CL.Patternfile_Count++;
				break;
			case 'G':
				CL.EXONGAP=atoi(optarg);
				break;
			case 'o':
				CL.JUNCTIONFILE=optarg;
				break;
			case 'n':
				MAX_MISMATCHES=atoi(optarg);
				break;
			/*case 'j':
				DUMP_ALL_JUNC=TRUE;
				break;
			case'b':
				MARKEX=TRUE;
				break;
			case 'm':
				MAPMODE=TRUE;PROCESSMODE=FALSE;
				break;
			case 'M':
				MAX_HITS_TO_STORE=atoi(optarg);
				break;
			case 'p':
				MAPMODE=FALSE;PROCESSMODE=TRUE;
				break;
			case 'N':
				COUNT=atoi(optarg);
				break;
			case 'r':
				USEREFGENE=TRUE;
				break;
			case 'w':
				WIGGLEFILE=optarg;
				break;
			case 'O':
				WRITE_SPLITREAD=TRUE;
				if (optarg) MAPFILE=optarg;
				break;*/
			case 'g':
				Name=optarg;Last_Dash=0;Genome_Name=optarg;
				for(;Name[0]!=0;Name++)
				{
					if (Name[0]=='/') 
					{
						Last_Dash++;Genome_Name=Name;
					}
				}

				Ind.REVBWTINDEX = (char*)Source;
				if(Last_Dash) Last_Dash=Genome_Name-optarg+1; else Genome_Name--;
				strncpy(Ind.REVBWTINDEX,optarg,Last_Dash);
				Ind.REVBWTINDEX[Last_Dash+0]='r';Ind.REVBWTINDEX[Last_Dash+1]='e';Ind.REVBWTINDEX[Last_Dash+2]='v';
				strcpy(Ind.REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(Ind.REVBWTINDEX+Last_Dash+3,".bwt"); 

				Ind.BWTFILE=Ind.REVBWTINDEX+500;
				strncpy(Ind.BWTFILE,optarg,Last_Dash);
				strcpy(Ind.BWTFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.BWTFILE+Last_Dash,".bwt"); 


				Ind.REVOCCFILE = Ind.BWTFILE+500;
				strncpy(Ind.REVOCCFILE,optarg,Last_Dash);
				Ind.REVOCCFILE[Last_Dash+0]='r';Ind.REVOCCFILE[Last_Dash+1]='e';Ind.REVOCCFILE[Last_Dash+2]='v';
				strcpy(Ind.REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(Ind.REVOCCFILE+Last_Dash+3,".fmv"); 


				Ind.OCCFILE=Ind.REVOCCFILE+500;			
				strncpy(Ind.OCCFILE,optarg,Last_Dash);
				strcpy(Ind.OCCFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.OCCFILE+Last_Dash,".fmv"); 

				Ind.SAFILE=Ind.OCCFILE+500;			
				strncpy(Ind.SAFILE,optarg,Last_Dash);
				strcpy(Ind.SAFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.SAFILE+Last_Dash,".sa");

				Ind.REVSAFILE = Ind.SAFILE+500;
				strncpy(Ind.REVSAFILE,optarg,Last_Dash);
				Ind.REVSAFILE[Last_Dash+0]='r';Ind.REVSAFILE[Last_Dash+1]='e';Ind.REVSAFILE[Last_Dash+2]='v';
				strcpy(Ind.REVSAFILE+Last_Dash+3,Genome_Name+1);
				strcat(Ind.REVSAFILE+Last_Dash+3,".sa"); 

				Ind.BINFILE=Ind.REVSAFILE+500;			
				strncpy(Ind.BINFILE,optarg,Last_Dash);
				strcpy(Ind.BINFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.BINFILE+Last_Dash,".pac");

				Ind.LOCATIONFILE=Ind.BINFILE+500;			
				strncpy(Ind.LOCATIONFILE,optarg,Last_Dash);
				strcpy(Ind.LOCATIONFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.LOCATIONFILE+Last_Dash,".ann.location");

                                Ind.BLKFILE = Ind.LOCATIONFILE+500;
                                strncpy(Ind.BLKFILE,optarg,Last_Dash);
                                strcpy(Ind.BLKFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.BLKFILE+Last_Dash,".blk.");
				sprintf(Ind.BLKFILE+strlen(Ind.BLKFILE),"%d",RQFACTOR);

                                Ind.INDFILE = Ind.BLKFILE+500;
                                strncpy(Ind.INDFILE,optarg,Last_Dash);
                                strcpy(Ind.INDFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.INDFILE+Last_Dash,".ind.");
				sprintf(Ind.INDFILE+strlen(Ind.INDFILE),"%d",RQFACTOR);

                                Ind.RANGEFILE = Ind.INDFILE+500;
                                strncpy(Ind.RANGEFILE,optarg,Last_Dash);
                                strcpy(Ind.RANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.RANGEFILE+Last_Dash,".range");

                                Ind.SORTEDRANGEFILE = Ind.RANGEFILE+500;
                                strncpy(Ind.SORTEDRANGEFILE,optarg,Last_Dash);
                                strcpy(Ind.SORTEDRANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.SORTEDRANGEFILE+Last_Dash,".sort");

                                Ind.INFOFILE = Ind.SORTEDRANGEFILE+500;
                                strncpy(Ind.INFOFILE,optarg,Last_Dash);
                                strcpy(Ind.INFOFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.INFOFILE+Last_Dash,".info");
				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
	CL.Patternfile_Count--;
	/*if (MARKEX) {Mark_Exons();exit(0);}*/
}

void Load_FM_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool)
{
	fprintf (stderr,"Loading index %s\n",Genome_Files.BWTFILE);
	fwfmi= Load_Indexes(Genome_Files.BWTFILE,Genome_Files.OCCFILE,Genome_Files.SAFILE,mmPool);
	fprintf (stderr,"Loading index %s\n",Genome_Files.REVBWTINDEX);
	revfmi= Load_Indexes(Genome_Files.REVBWTINDEX,Genome_Files.REVOCCFILE,Genome_Files.REVSAFILE,mmPool);
	fwfmi->saInterval=revfmi->saInterval;
	fprintf (stderr,"Loading Location %s\n",Genome_Files.LOCATIONFILE);
	int Genome_Count= Load_Location(Genome_Files.LOCATIONFILE,Annotations,NULL);
	/*fprintf (stderr,"Loading packed genome %s\n",Genome_Files.BINFILE);
	loadPac(Genome_Files.BINFILE);*/
	fprintf (stderr,"Loading index %s\n",Genome_Files.INDFILE);
	fprintf(stderr,"Done...\n");
}


void InitX(BWT *revfmi,unsigned & SOURCELENGTH,gzFile & Input_File,gzFile & Mate_File,FILETYPE & File_Info,Parameters & CL,Index_Info & Genome_Files,int & INIT_MIS_SCAN)
{
	EXONGAP=CL.EXONGAP;
	SOURCELENGTH = revfmi->textLength;
	Char_To_Code['N']=0;Char_To_Code['n']=0;
	Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;
	Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;
	Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;
	Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;
	Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;
	Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt

	Open_Files(Input_File,Mate_File,CL);
	Detect_Input(File_Info,Input_File,Mate_File);

	CONVERSION_FACTOR=revfmi->textLength-File_Info.STRINGLENGTH;
	Conversion_Factor=revfmi->textLength-File_Info.STRINGLENGTH;
	
	if(INIT_MIS_SCAN==-1)
	{
		if (File_Info.STRINGLENGTH<=80)
		{
			INIT_MIS_SCAN=1;
		}
		else if (File_Info.STRINGLENGTH<=100)
		{
			INIT_MIS_SCAN=2;
		}
		else
		{
			INIT_MIS_SCAN=3;
		}
	}
	cout <<"Scanning @ " <<INIT_MIS_SCAN<<"Mismatches\n";
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Convert Reverse
 *  Description:  Convert Read in binary form to its binary and nucleotide rev. comp.
 * =====================================================================================
 */

void Convert_Reverse(char* Read,char * RC_Read,char* RC_Bin,int StringLength)
{
	for (unsigned i=0;i<=StringLength-1;i++)
	{
		RC_Bin[StringLength-1-i]=3-Read[i];
		RC_Read[StringLength-1-i]="ACGT"[3-Read[i]];
	}
}

bool  Progress_Bar(Parameters & CL,unsigned & Number_of_Tags,unsigned & Progress,unsigned & Tag_Count,FILETYPE & File_Info)
{
	Tag_Count++;
	Progress++;
	if (CL.MAX_TAGS_TO_PROCESS && Tag_Count >= CL.MAX_TAGS_TO_PROCESS) return false; 

	if (Progress==Number_of_Tags) 
	{
		if (CL.MAX_TAGS_TO_PROCESS)
		{
			Number_of_Tags=(CL.MAX_TAGS_TO_PROCESS)/20;
			Progress=0;
			Show_Progress(Tag_Count*100/CL.MAX_TAGS_TO_PROCESS);
		}
		else
		{
			off64_t Current_Pos=ftello64(File_Info.Org_File);
			unsigned Average_Length=Current_Pos/Tag_Count+1;//+1 avoids divide by zero..
			Number_of_Tags=(File_Info.File_Size/Average_Length)/20;
			Progress=0;
			off64_t Perc= Current_Pos*100/File_Info.File_Size;
			Show_Progress(Perc);
		}
	}
	return true;
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

int Map_Read(MEMX & MF,MEMX & MC,int MAX_MISMATCHES, LEN & L,BWT* fwfmi,BWT* revfmi,int Next_Mis,int Max_Hits)
{
	if(MAX_MISMATCHES < Next_Mis) return -1;
	assert(Next_Mis >=0);assert(MAX_MISMATCHES >= Next_Mis);assert (Next_Mis <= 5);
	int In_Mis=0,Hits=0;MF.Hits=0;//MC.Hits=0;
	if (Next_Mis == 0) goto Zero; else if (Next_Mis ==1) goto One;else if (Next_Mis ==2) goto Two;else if (Next_Mis ==3) goto Three;else if (Next_Mis ==4) goto Four; else goto Five;
	assert(Next_Mis||MF.Hit_Array_Ptr);
Zero:
	assert(MF.Hit_Array_Ptr==0);
	Hits+=Zero_Mismatch(MF.Current_Tag,L,revfmi,MF);
	if(Hits<Max_Hits)
		Hits+=Zero_Mismatch(MC.Current_Tag,L,revfmi,MC);
One:
	if (!Hits && MAX_MISMATCHES >0)
	{
		In_Mis=1;
		Hits+=One_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		if(Hits<Max_Hits)
			Hits+=One_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Two:
	if (!Hits && MAX_MISMATCHES >1)
	{
		In_Mis=2;
		Hits+=Two_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		if(Hits<Max_Hits)
			Hits+=Two_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Three:
	if (!Hits && MAX_MISMATCHES >2)
	{
		In_Mis=3;
		Hits+=Three_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		if(Hits<Max_Hits)
			Hits+=Three_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Four:
	if (!Hits && MAX_MISMATCHES >3)
	{
		In_Mis=4;
		Hits+=Four_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		if(Hits<Max_Hits)
			Hits+=Four_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Five:
	if (!Hits && MAX_MISMATCHES >4)
	{
		In_Mis=5;
		Hits+=Five_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		if(Hits<Max_Hits)
			Hits+=Five_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}

	MF.Hit_Array[MF.Hit_Array_Ptr].Start=0;//MC.Hit_Array[MC.Hit_Array_Ptr].Start=0;//tag sentinels to sa lists..
	//MF.Hit_Array_Ptr++;//MC.Hit_Array_Ptr++;//Setup for suboptimal hits..
	assert(In_Mis <= MAX_MISMATCHES);
	//Top=Hits;
	return (Hits ? In_Mis : -1) ;
}

void Open_Outputs(ofstream & SAM,string filename)
{
	SAM.open(filename.c_str());
	if(SAM.is_open())
	{
		return;
	}
	else
	{
		cout <<"File open error\n";
		exit(100);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------

bool Call_Junc(char* Junc,unsigned Pos,int StringLength,ofstream & Out_File,string Des)
{
	//assert (L.c_str());// && DesS.c_str());
	char Dummy[5000],Right[500],Left[500],Chr[20],Strand[3];//,Read[500],SignM;
	unsigned RightL,LeftL,RightR,LeftR;

	/*char *token = strtok (Junc, ",");
	sscanf (token, "%s", Left);
	token = strtok (NULL,",");
	if(token)
	{
		    sscanf (token, "%s", Right);
	}
	else {Right[0]=0;}*/

	//char *token;
	char *token= strchr (Junc, ',');
	if(!token) 
		return 0;
	*token=0;
	sscanf (Junc, "%s", Left);
	sscanf (token+1, "%s", Right);

	if(!Right[0] || 'E'==Right[0]|| '<'==Right[0])
	{
		return 0;
	}
	else
	{
		sscanf(Left,"%[^:]%*c%u%*c%u%s",Chr,&LeftL,&LeftR,Strand);
		sscanf(Right,"%*[^:]%*c%u%*c%u%s",&RightL,&RightR,Strand);
		char Sign=Strand[1];
		assert(LeftL>0 && LeftR >0 && RightL >0 && RightR>0 && (Sign=='+' || Sign=='-'));

		int Left_Gap=LeftR-(LeftL+Pos);
		assert(Left_Gap >0);
		if ((Left_Gap<RESIDUE) || (Left_Gap>StringLength-RESIDUE))
		{
			return 0;
		}
		else
		{

			bool Skip_Hit=false;
			if(!Skip_Hit)
			{
				int Left_Cord=LeftR-3;int Right_Cord=RightL+3;
				int Gap=Right_Cord-Left_Cord-3;
				assert(Gap >0);
				Out_File \
					<<Chr<<"\t"<<Left_Cord<<"\t"<<Right_Cord<<"\t" \
					<<"JUNCXXX\t1000\t+\t" \
					<<Left_Cord<<"\t"<<Right_Cord-1 \
					<<"\t255,0,\t2\t3,3\t0," \
					<<Gap<<endl;
				/*cout \
					<< Des <<"\t" \
					<<Left<<"\t"<<Right<<"\t"<<Pos<<endl;*/
				return 1;
			}
			else return 0;
		}
	}
	return true;
}
