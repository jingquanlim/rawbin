#ifndef COMMON_H
#define COMMON_H
#include "const.h"
#define SOLEXA_READS 0
#define ENT_LIM 50 
#define SOLID_READS  1
#define SINGLE_END 1
#define PAIRED_END 2
#define MAX_MISMATCHES_BOUND 16
#define REVERSE 11
#define FORWARD 10
#define SAINTERVAL 8
#define BRANCHTHRESHOLD 0//80 //30 //Threshold at which to check the BWT instead of branching
#define PAIREND 2
#define INDELMARK (MAXTAG-1)
#define INSERTMARK 63
#define DELETEMARK 70
#define BATFILEMODE 	1	//pairing from already mapped batman file..
#define NORMALFILEMODE  2
#include <limits.h>
#include "zlib.h"
//`#include <map>
#include "Hash.h"
extern "C" 
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}
int const MAX_CIGLEN=50;

//{-----------------------------  STRUCTS  -------------------------------------------------/
struct Ann_Info
{
	unsigned Size;
	unsigned Cumulative_Size;
	int ID;
	char* Name;
};

/*struct Offset_Record
{
	char Genome[40];
	unsigned Offset;
};*/

struct Offset_Record
{
	char Genome[40];
	char GenomeM[40];
	unsigned Offset;
	/*FILE* Out_File;
	FILE* Out_FileM;
	FILE* Unmapped;
	FILE* UnmappedX;
	FILE* UnmappedM;
	FILE* UnmappedXM;
	FILE* Ref_File;*/
	Hash*  Junc_Hash;
};

struct OUTPUT//info for writing output..
{
	char SAM;
	char PLUSSTRAND;
	char MaxHits;
	char Offset;
	char* Buffer;
	unsigned* Location_Array;
	Offset_Record* Genome_Offsets;
	int Genome_Count;
	int Length_Array[3];
	char FILETYPE;
};

struct Header
{
	char ID[3] ;
	unsigned MAXHITS;
	char FILETYPE;
	char HITMODE;
	char IGNOREHEAD;
	char Index_Count;
	int Tag_Length;
	char Print_Desc;//print the tag desc?
}__attribute__((__packed__));

struct In_File
{
	int MAXHITS;
	int STRINGLENGTH;
	char FILETYPE;
	char MAX_MISMATCHES;
	int Stat_Size;
	int TAG_COPY_LEN;
	char ROLLOVER;
	char SCANBOTH;
	char PRINT_DESC;
	char LOADREVERSEONLY;
	char NORMAL_TAGS;
	int Length_Array[3];
	int HEAD_LENGTH,TAIL_LENGTH;
	char *Positive_Head,*Positive_Tail;
	char Tag_Copy[MAXTAG+1];
	off64_t File_Size;
};


// LEN stores the information related to stringlength...
struct LEN
{
	int STRINGLENGTH_ORG;//Real length of string
	int STRINGLENGTH;//effective length of string = Forced length-ignore head
	int STRINGLENGTHl;
	int STRINGLENGTHr; 
	int LH; 
	int RH; 
	int LHQL; 
	int LHQR; 
	int RHQL; 
	int RHQR; 
	int LHl; 
	int RHl; 
	int LHQLl; 
	int LHQRl; 
	int RHQLl; 
	int RHQRl; 
	int LHr; 
	int RHr; 
	int LHQLr; 
	int LHQRr; 
	int RHQLr; 
	int RHQRr;
	int IGNOREHEAD;
};

struct READ
{
	char Tag_Copy[2*MAXTAG+2];
	char Description[MAXDES+1];
	char Tag[MAXTAG+1];
	char Quality[MAXTAG+1];
	char Plus[MAXTAG+1];
	char NLocations[MAXTAG+1];
	char Complement[MAXTAG+1];
	int NCount;
	unsigned char N[500];
};

struct SAMREAD
{
	char Chr[20];
	unsigned Loc;
	int MapQ; 
	char Cigar[MAX_CIGLEN+1];
	char SAM_Line[5000];
	bool Enable;
	int NM;
};

/*struct READ
{
	char Description[MAXDES];
	char Tag_Copy[MAXTAG];
	char Quality[MAXTAG];
	char Plus[MAXTAG];
	int NCount;//Number of N's
	char N[MAXTAG];
	char NLocations[MAXTAG];
	char Tag_Number;//Head =1, Tail =2
	unsigned Read_Number;
};*/

struct FMFILES
{
	char* DISCORDANTFILE;
	char* PATTERNFILE; 
	char* PATTERNFILE1;
	char* HITSFILE;//=HITSFILE_DEF;//file to output hits...
	char* BWTFILE ; 
	char* OCCFILE ;
	char* REVBWTINDEX;
	char* REVOCCFILE;
	char* REVSAFILE;
	char* SAFILE;
	char* BINFILE;
	char* PACFILE;
	char* LOCATIONFILE ; 
	char* INPUTFILE;
	char* OUTPUTFILE;
	char* INDFILE;
	char* RANGEFILE;
	char* NLOCATIONFILE;
	char* BLKFILE;
	char* SINGLEFILE;

};

struct BATREAD
{
	int StringLength;
	int IGNOREHEAD;
	int NCount;
	char Forward[MAXTAG];
	char Complement[MAXTAG]; 
};

struct INFILE
{
	FILE* Input_File;
	off64_t File_Size;
	char FILETYPE;
	char SOLID;
	int TAG_COPY_LEN;
	char TAB_SEPERATED;
	char *Buffer;
	int PAIR_LENGTH_RIGHT;
	int PAIR_LENGTH_LEFT;
};

struct FILELIST
{
	FILE* Head;
	FILE* Tail;
};

struct BATPARAMETERS
{
	unsigned MAXHITS;
	char* PATTERNFILE;
	char* PATTERNFILE1;
	char* MISFILE1;
	char* MISFILE2;
	char* UNMAPPED_FILE;
	char* MULTI_FILE;
	char ONEFMINDEX;
	char UNMAPPED;
	char MAX_MISMATCHES;
	char PAIRING_MODE;
	unsigned MAX_TAGS_TO_PROCESS;
	int NTHREADS;
	int IGNOREHEAD; 
	int FORCELENGTH;
	int INSERTSIZE;
	unsigned OFFSET;
	int PLUSSTRAND;
	int SW_FLANKSIZE;//Length of string to do S/W... A portion of this size will be taken from the either side of insert size + good hit.
	int FLANKSIZE;
	char LOADREVERSEONLY;
	char MISHITS;
	char USELOCATION; 
	int ORIENTATIONS;//how many strand combinations to try..
	int Patternfile_Count;//0=single end,1= paired end.. 
	int Misfile_Count;
	char ROLLOVER;
	char SCANBOTH;
	char SMITH_WATERMAN;
	char VERIFY;
	char LOG;
	char SOLIDMAP;
	int STD;
	char FORCESOLID;
};

struct FastaSeq 
{
      int   len; /* the actual string length of seq */
      char* seq; /* the sequence itself */
};

typedef struct gz_stream {
    z_stream stream;
    int      z_err;   /* error code for last stream operation */
    int      z_eof;   /* set if end of input file */
    FILE     *file;   /* .gz file */
    Byte     *inbuf;  /* input buffer */
    Byte     *outbuf; /* output buffer */
    uLong    crc;     /* crc32 of uncompressed data */
    char     *msg;    /* error message */
    char     *path;   /* path name for debugging only */
    int      transparent; /* 1 if input file is not a .gz file */
    char     mode;    /* 'w' or 'r' */
    z_off_t  start;   /* start of compressed data in file (header skipped) */
    z_off_t  in;      /* bytes into deflate or inflate */
    z_off_t  out;     /* bytes out of deflate or inflate */
    int      back;    /* one character push-back */
    int      last;    /* true if push-back is last character */
} gz_stream;


struct HEADER
{
	char ID[3] ;
	unsigned MAXHITS;
	char FILETYPE;
	char HITMODE;
	char IGNOREHEAD;
	char Index_Count;
	int Tag_Length;
	char Print_Desc;//print the tag desc?
}__attribute__((__packed__));

struct SARange
{

	unsigned long Start;
	unsigned long End;
	int Level;//at what relative node are we?
	char Mismatches;//number of mismatches at this level
	unsigned char Skip;
	//int Tag;//Tag number
	char FMIndex;
	char Strand;
	unsigned Mismatch_Char;//2|2|...
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
        unsigned Mismatch_PosX;
};

struct Range
{
	unsigned Start;
	unsigned End;
	int Label;//Final Label of the range...
};

struct Mismatches_Record
{

	int 	Gap;
	unsigned Mismatch_PosX;//int?
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
	unsigned Mismatch_Char;//2|2|...


}__attribute__((__packed__));

struct Mismatches_Record_GIS
{

	unsigned Mismatch_Char;//2|2|...
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...

}__attribute__((__packed__));


struct Output_Record
{
	unsigned Tag;
	unsigned  Start;
	char Index;
	unsigned char Skip;
	char Mismatches;
	int Gap;
}__attribute__((__packed__));

struct Branches
{
	long Is_Branch [4];
};

struct MEMLOOK
{
	int Lookupsize;
	unsigned* Forward_Start_LookupX;
	unsigned* Forward_End_LookupX;
	unsigned* Backward_Start_LookupX;
	unsigned* Backward_End_LookupX;
};

/*struct MEMX
{
	SARange Branch_Ranges[4];

	unsigned* Forward_Start_LookupX;
	unsigned* Forward_End_LookupX;
	unsigned* Backward_Start_LookupX;
	unsigned* Backward_End_LookupX;
	unsigned ARRAY_BOUND;
	unsigned END_BOUND;
	unsigned Hits;
	unsigned short Stats[7];

	char Guessed;
	char FMIndex;
	char Strand;
	char Stat_Size;
	char Larger_Than_Ten;

	char* Write_Buffer;
	char* Current_Tag;

	int Hit_Array_Ptr;
	int Lookupsize;
	int Last_Mismatch_Written;
	SARange* Hit_Array;
	READ Read;
	OUTPUT Output;

	int Left_Mishits_Pointer;
	int Right_Mishits_Pointer;
	int Possible_20_Pointer;
	int Possible_02_Pointer;
	int Mismatches_Forward_Pointer;
	int Mismatches_Backward_Pointer;
	int Two_Mismatches_At_End_Pointer;
	int Two_Mismatches_At_End_Forward_Pointer;
	int Possible_03_Pointer;
	int Possible_30_Pointer;
	int Possible_04_Pointer,Possible_40_Pointer,Possible_50_Pointer;
	int Possible_05_Pointer;
	int Mismatches_Forward_Pointer_Last4;
	int Left_Mishits_Pointer_1;
	int Mismatches_Forward_Pointer_Last5;

	SARange* BMHStack;
	SARange* FSHStack;
	SARange* FSHStackX0X;
	SARange* FSSStack;
	SARange* FSSStackX;
	SARange* BMStack;
	SARange* BMStackX;
	SARange* BMStack_X11;
	SARange* BMStack_X11H;
	SARange* PSBStack;

	SARange* Exact_Match_Forward;
	SARange* Exact_Match_Backward;
	SARange* Left_Mishits;
	SARange* Right_Mishits;
	SARange* Mismatches_Backward;
	SARange* Mismatches_Forward;
	SARange* Two_Mismatches_At_End_Forward;
	SARange* Two_Mismatches_At_End;
	SARange* Possible_20;
	SARange* Possible_02;
};*/
struct MEMX
{
	SARange Branch_Ranges[4];

	unsigned* Forward_Start_LookupX;
	unsigned* Forward_End_LookupX;
	unsigned* Backward_Start_LookupX;
	unsigned* Backward_End_LookupX;
	unsigned ARRAY_BOUND;
	unsigned END_BOUND;
	unsigned Hits;
	unsigned short Stats[7];

	char Guessed;
	char FMIndex;
	char Strand;
	char Stat_Size;
	char Larger_Than_Ten;
	char Extend;

	char* Write_Buffer;
	char* Current_Tag;

	int Hit_Array_Ptr;
	int Lookupsize;
	int Last_Mismatch_Written;
	SARange* Hit_Array;
	READ Read;
	OUTPUT Output;
	LEN L;

	int Left_Mishits_Pointer;
	int Right_Mishits_Pointer;
	int Possible_20_Pointer;
	int Possible_02_Pointer;
	int Mismatches_Forward_Pointer;
	int Mismatches_Backward_Pointer;
	int Two_Mismatches_At_End_Pointer;
	int Two_Mismatches_At_End_Forward_Pointer;
	int Possible_03_Pointer;
	int Possible_30_Pointer;
	int Possible_04_Pointer,Possible_40_Pointer,Possible_50_Pointer;
	int Possible_05_Pointer;
	int Mismatches_Forward_Pointer_Last4;
	int Left_Mishits_Pointer_1;
	int Mismatches_Forward_Pointer_Last5;

	int Least_Mis;
	int Best_Quality;

	SARange* BMHStack;
	SARange* FSHStack;
	SARange* FSHStackX0X;
	SARange* FSSStack;
	SARange* FSSStackX;
	SARange* BMStack;
	SARange* BMStackX;
	SARange* BMStack_X11;
	SARange* BMStack_X11H;
	SARange* PSBStack;

	SARange* Exact_Match_Forward;
	SARange* Exact_Match_Backward;
	SARange* Left_Mishits;
	SARange* Right_Mishits;
	SARange* Mismatches_Backward;
	SARange* Mismatches_Forward;
	SARange* Two_Mismatches_At_End_Forward;
	SARange* Two_Mismatches_At_End;
	SARange* Possible_20;
	SARange* Possible_02;
};

struct GUESS
{
	MEMX Guessed;
	char* Guessed_Read;
	MEMX Guess_Complement;
	char* Guessed_ReadC;
};
struct MEM
{
	unsigned* Forward_Start_LookupX;
	unsigned* Forward_End_LookupX;
	unsigned* Backward_Start_LookupX;
	unsigned* Backward_End_LookupX;
	unsigned ARRAY_BOUND;
	unsigned END_BOUND;

	char* Write_Buffer;
	int Lookupsize;

	SARange* BMHStack;
	SARange* FSHStack;
	SARange* FSHStackX0X;
	SARange* FSSStack;
	SARange* FSSStackX;
	SARange* BMStack;
	SARange* BMStackX;
	SARange* BMStack_X11;
	SARange* BMStack_X11H;
	SARange* PSBStack;

	SARange* Exact_Match_ForwardF;
	SARange* Exact_Match_BackwardF;
	SARange* Left_MishitsF;
	SARange* Right_MishitsF;
	SARange* Mismatches_BackwardF;
	SARange* Mismatches_ForwardF;
	SARange* Two_Mismatches_At_End_ForwardF;
	SARange* Two_Mismatches_At_EndF;
	SARange* Possible_20F;
	SARange* Possible_02F;

	SARange* Exact_Match_ForwardC;
	SARange* Exact_Match_BackwardC;
	SARange* Left_MishitsC;
	SARange* Right_MishitsC;
	SARange* Mismatches_BackwardC;
	SARange* Mismatches_ForwardC;
	SARange* Two_Mismatches_At_End_ForwardC;
	SARange* Two_Mismatches_At_EndC;
	SARange* Possible_20C;
	SARange* Possible_02C;

};


struct Thread_Arg
{
	BWT* fwfmi;
	BWT* revfmi;
	FILE* Output_File;
	LEN  L; 
	char ONEFMINDEX;
	MEMLOOK MLook;
	OUTPUT Output;
	BATPARAMETERS Bat_Para;
	INFILE In;
	int ThreadID;
};


struct Threading
{
	pthread_t Thread;
	unsigned r;
	void *ret;
	Thread_Arg Arg;
};
//}-----------------------------  STRUCTS  -------------------------------------------------/*
#include <signal.h>
#define Breakpoint raise(SIGINT)

#endif
