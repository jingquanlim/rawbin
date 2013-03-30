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
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
//#include <dvec.h>
#include "bfix.h"
#include <getopt.h>
#include "zlib.h"
#include "assert.h"
#include "const.h"
//#include "rqindex.h"
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

using namespace std;
//{-----------------------------  DEFINES  -------------------------------------------------/

typedef map <unsigned,unsigned> HASH_PAIR; 
typedef map <unsigned,unsigned char> JUNC_HASH;
#define EXTEND STRINGLENGTH/2//8 
const int MAX_JUNC_LIST=500000;
const int EXTENDISLAND=10;
const int OVERLAP=10;
const int MAX_HITS_ALLOWED=2;//maximum pairs to output..
const int GAPLIMIT=7;//maximum pairs to output..
const int COVERAGE_THRESH=0;//5
const int RESIDUE=6;
const int MAXMULTIGAP=50000;
const int HIGHREPTHRES = 30;
const int EX_MIN=4;//10;//How close the junctions are..
const int REFMARK=10000;
const int MAX_REF_LINE=10000;
const int MAX_X=500;//for mouse 312
const int MAX_ONE_SIDE_HITS=5;
const int MAXSTRINGLEN=255;
const int SKIP_INIT_HIT=1;//If the full read can be mapped, dont scan for junctions ...
const int MAX_JUNC_IN_READ=100;//If the full read can be mapped, dont scan for junctions ...
enum { FULLREAD,SPLITREAD,UNIQUE,UNIQUEX,NON_UNIQUE,UNMARKED,UNEXTENDED};

#define MIS_IN_UNMAPPED 0
#define MIS_IN_RES 0
#define N_a	0
#define N_c	1
#define N_g	2
#define N_t	3
#define S_non	0
#define S_dummy	500
#define S_ag	1
#define S_gt	2
#define S_ac	3
#define S_gc	4
#define S_at	5
#define S_plus  255
#define S_gt_ag 256
#define S_gc_ag 257
#define S_at_ac 258

#define S_ct	6
#define S_ct_ac 259
#define S_ct_gc 260
#define S_gt_at 261
#define S_minus  262
#define POS_JUNC 	256
#define MINUS_JUNC 	259

#define INIT 100

#define LEFT	0
#define LEFTC	1
#define RIGHT	2
#define RIGHTC	3

#define FW	1
#define BW	0
#define MAX_REFGENE_OVERHANG 20	//maximum overhang allowed to determine for overhang over a putative splice.

#define INREFGENE	1
#define DONOR		2
#define ACCEPTOR	4

#define PLUS	0
#define MINUS	1
#define DELETED 2
#define EXON	4
#define COVTYPE unsigned char //what type can the maximum coverage be?
#define COVTYPEMAX UCHAR_MAX //what type can the maximum coverage be?
//#define COVTYPE unsigned short //what type can the maximum coverage be?
//#define COVTYPEMAX USHRT_MAX //what type can the maximum coverage be?

#define TAB	1
#define FQ	2
#define FA	3
#define SAINTERVAL 8
#define SAGAP_CUTOFF 5
#define SAGAP_CUTOFF_INDEX 2
#define TWOFILE	4
#define GISMODE 100
#define MAXTAG 256
#define MAXDES 500
#define PAIR_END_SEPERATOR '\t'
#define MAX_MISMATCHES_BOUND 2
#define BRANCHTHRESHOLD 80 //30 //Threshold at which to check the BWT instead of branching
#define REVFMI	1 
#define FWFMI	0 
#define BUFSIZE 100000
//#define MAXCOUNT 3 //30000  //Maximum number of pairs to enum...
//#define MAX_HITS_TO_STORE 200//later resolve this and MAXCOUNT conflict..
//#define MAXCOUNT MAX_HITS_TO_STORE //30000  //Maximum number of pairs to enum...
#define HEAD 0
#define TAIL 1

#define _lrotr(x, n) ((((unsigned long)(x)) >> ((int) ((n) & 31))) | (((unsigned long)(x)) << ((int) ((-(n)) & 31))))
#define _lrotl(x, n) ((((unsigned long)(x)) << ((int) ((n) & 31))) | (((unsigned long)(x)) >> ((int) ((-(n)) & 31))))
#define ror32(x, n) _lrotr(x, n)
#define rol32(x, n) _lrotl(x, n)
#define bswap32(x) (rol32((unsigned long)(x), 8) & 0x00ff00ff | ror32 ((unsigned long)(x), 8) & 0xff00ff00)


#ifdef BAT64
	#define UINT uint64_t	
	#define UINTSIZE 64
	#define MULTISTRINGLENTH 28
	//#define BSWAP __builtin_bswap64
	#define BSWAP bswap32
	#define SEARCH_MASK 0xFFFFFFFFFFFFFF00
	#define SEARCH_MASKB 0x00FFFFFFFFFFFFFF
	#define SEARCH_ISOLATE_SHIFT 62
	#define BSF __builtin_clzl
	#define BSR __builtin_ctzl
	#define INTLENGTH 64 
	#define BYTES_IN_UNIT  UINTSIZE/8
#else
	#define UINT uint32_t
	#define UINTSIZE 32 
	#define MULTISTRINGLENTH 12
	//#define BSWAP __builtin_bswap32
	#define BSWAP bswap32
	#define SEARCH_MASK 0xFFFFFF00
	#define SEARCH_MASKB 0x00FFFFFF
	#define SEARCH_ISOLATE_SHIFT 30
	#define BSF __builtin_clz
	#define BSR __builtin_ctz
	#define INTLENGTH 32
	#define BYTES_IN_UNIT  UINTSIZE/8
#endif
//}-----------------------------  DEFINES  -------------------------------------------------/

//{-----------------------------  STRUCTS  -------------------------------------------------/

struct OP//Ordered pairs..
{
	unsigned x;
	unsigned y;
	//int Motif;
};

struct OPX//Ordered pairs..
{
	unsigned x;
	unsigned y;
	int Motif;
};

struct JStat//Junction stats..
{
	int Count;
	int Junc_Type;
};

struct OP_Cmp
{
	bool operator()( OP OP1, OP OP2)
	{
		if (OP1.x == OP2.x)
		{
			return OP1.y< OP2.y;
		}
		else return OP1.x < OP2.x;
	} 
};

struct Out_Record
{
	unsigned ID;//ID for read number..
	char Strlen;
	int HLength;
	unsigned HLocation;
	int TLength;
	unsigned TLocation;
	int Junc_Count;// # of junctions in this read..
	char Type;
};

struct Unmapped_Record
{
	int Strlen;
	int Length;
	unsigned Location;
//next comes string...
};

struct State
{
	char State;
	unsigned Location;
	int Score;
};

struct Xon
{
	unsigned Start;
	unsigned End;
	int Start_Type;
	int End_Type;
};

struct Island
{
	unsigned Start;
	unsigned End;
	Island* Next;
	char Status;
};

struct Offset_Record
{
	char Genome[40];
	char GenomeM[40];
	unsigned Offset;
	FILE* Out_File;
	FILE* Out_FileM;
	FILE* Unmapped;
	FILE* UnmappedX;
	FILE* UnmappedM;
	FILE* UnmappedXM;
	FILE* Ref_File;
};

struct PAIR
{
	unsigned Head;
	unsigned Tail;
	int HLevel;
	int TLevel;
};

struct TAG_INFO
{
	unsigned SA_Start;//start of the Sa range
	unsigned Gap;//length of SA range
	unsigned Block_Start;//Start of block info..
	unsigned Index;//Index to SA_index
	unsigned First;//first location of hit
	unsigned Last;//last location of hit
	unsigned Field_Length;
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

struct SA
{

	unsigned Start;
	unsigned End;
	unsigned Start_Location;// exact location of the first occurance...
	unsigned End_Location;

};

struct SARANGE
{
	unsigned Start;
	unsigned End;
	int Level;//at what relative node are we?
	char Mismatches;//number of mismatches at this level
	unsigned char Skip;
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
	unsigned Mismatch_Char;//2|2|...
	unsigned Conversion_Factor;
	
};

struct ORIENTATION
{

	char Strand;
	char* String;
	SARANGE TagB;//RQfactor tag..
	SARANGE TagBR;//RQfactor tag..
	SARANGE TagBH;//RQfactor tag..
	SARANGE TagF;//Halfway tag..
	SARANGE TagFR;//Halfway tag..
	SARANGE TagFH;//Halfway tag..

};

struct RANGEINDEX
{
	FILE* Index;
	FILE* Blocks;
	SA *SA_Index;
	int *SA_Blocks;
	unsigned Hits;
	char COMPRESS;
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
//}-----------------------------  STRUCTS  -------------------------------------------------/*

//{-----------------------------  Classes  -------------------------------------------------/

typedef map <OP,JStat>::iterator map_it; 
class Hash
{

	public:
	map <OP,JStat,OP_Cmp> Junctions;
	map_it Junc_I,JJ;
	map_it Last;

	void Insert (OP Location,int Paring);
	void Delete (OP Location);
	bool Init_Iterate(OP & Location,JStat & Data);
	bool Iterate(OP & Location,JStat & Data);
	map_it Begin();


};


map_it Hash::Begin ()
{
	return Junctions.begin();
}

void Hash::Insert (OP Location,int Paring)
{
	JStat T;
	assert(Paring >=255 && Paring <=262); assert(Location.x <= Location.y);
	Junc_I=Junctions.find(Location);
	if (Junc_I == Junctions.end()) //New entry..
	{
		T.Count=1;T.Junc_Type=Paring;
		Junctions[Location]=T;
	}
	else 
	{
		assert((Junc_I->second).Count>0 && (Junc_I->second).Junc_Type==Paring);
		(Junc_I->second).Count++;
	}
}


void Hash::Delete (OP Location)
{
	Junc_I=Junctions.find(Location);
	if (Junc_I != Junctions.end()) Junctions.erase(Junc_I);
}

bool Hash::Init_Iterate(OP & Location,JStat & Data)
{
	Junc_I=Junctions.begin();
	if (Junc_I == Junctions.end()) return false; 
	else 
	{
		Last=Junc_I;//save last pos...
		Location = Junc_I->first;
		Data = Junc_I->second;
		Junc_I++;
		return true;
	}
}

bool Hash::Iterate(OP & Location,JStat & Data)
{
	if(Junc_I == Junctions.end()) return false;
	Last=Junc_I;//save last pos...
	Location = Junc_I->first;
	Data = Junc_I->second;
	Junc_I++;
	return true;
}

//}-----------------------------  Classes  -------------------------------------------------/

//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

bool Scan(unsigned & Total_Hits,unsigned & Tags_Processed,int Off,unsigned Read_ID);
char* Get_Junc_Type(int Junc_Type);
void Load_RefGene(FILE* Ref_Handle,Hash & Ref_Gene,JUNC_HASH & Donor_AT,JUNC_HASH & Donor_GT,JUNC_HASH & Donor_GC,JUNC_HASH & Donor_CT,JUNC_HASH & Acc_AC,JUNC_HASH & Acc_AG,JUNC_HASH & Acc_AT,JUNC_HASH & Acc_GC);
void Init();
void Allocate_Memory();
void Open_Files();
int Split(char* String,char Sep, char* Fields[],int Max=0);
void Load_Range_Index(char* INDFILE, char* BLKFILE,RANGEINDEX & Range_Index);
void Parse_Command_line(int argc, char* argv[]);
void Detect_Input(int & TAG_COPY_LEN,int & STRINGLENGTH,char & NORMAL_TAGS, char & PAIRING_TYPE,char & FILETYPE,int & PAIR_LENGTH_RIGHT,int & PAIR_LENGTH_LEFT);
void Convert_To_REVSA(SARANGE & Tag, char* Current_Tag);
void Get_Head_Tail(SARANGE & Head, SARANGE & Tail,unsigned d,PAIR* Pairs,int & Pairs_Index);
void Search_Small_Gap(SARANGE & Head, SARANGE & Tail, unsigned d,PAIR* Pairs,int & Pairs_Index);
void Load_Info( TAG_INFO & Tag, SARANGE & Head);
void Print_BED(Hash & Junctions,Island* Island_List, COVTYPE* Coverage,char* chromosome,char Strand);
char Two_Mismatch_Search();
void Show_Progress(unsigned Percentage);
char Search_Forwards(SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Forward_Pointer,SARANGE* Two_Mismatches_At_End_Forward,int & Mismatches_Forward_Pointer,SARANGE* Mismatches_Forward);
char Search_Backwards(struct SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int & Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward);
void Search_Forwards_OneSA(SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Forward_Pointer,SARANGE* Two_Mismatches_At_End_Forward,int & Mismatches_Forward_Pointer,SARANGE* Mismatches_Forward);
void Branch_Detect_Backwards (const struct SARANGE Tag,BWT *fmi,int Start);
void Branch_Detect (const struct SARANGE Tag,BWT *fmi,int Start);
unsigned Get_Location(TAG_INFO & Tag, unsigned Offset);
void Pair_Reads(SARANGE *Head_Hits,SARANGE *Tail_Hits);
int Extend_Location_Forward(SARANGE & Tag,int Start,int StringLength,int Count);
int Extend_Location_Backward(SARANGE & Tag,int Start,int StringLength,int Count);
void Verbose();
void Pack_Text(char* Current_Tag);
unsigned Get_File_Size(FILE* File);
int Location_To_Genome(unsigned & Location);
int Load_Locations(Offset_Record* Genome_Offsets,unsigned Offsets[]);
void Get_SARange_Fast( char New_Char, SARANGE & Range,BWT *fmi);
void Coverage(Offset_Record Current_Genome);
void Print_LocationX(SARANGE & Tag,char Index);
void Write_Sam(char Sign,char* Chr,unsigned Loc1,unsigned Loc2,int G1,int G2,char Code,unsigned Read_ID);
int Juncs_In_Read(OPX* Junc_List,int & Junc_List_Ptr,int Head_Start,int Head_Length,int Tail_Start, int Tail_Length);
void Find_Split_Junctions(OPX* Junc_List,int Junc_List_Ptr,Hash & Junctions,unsigned Last_ID,char* Genome);

int Scan_Minus_Motif(unsigned x,unsigned y);
int Scan_Plus_Motif(unsigned x,unsigned y);
void Mark_Exons();
void Seek_Island_Junctions(Island* Island_List,map <unsigned,unsigned char> & Donor_AT,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Donor_GC,map <unsigned,unsigned char> & Acc_AG,map <unsigned,unsigned char> & Acc_AC);
void Seek_Island_JunctionsM(Island* Island_List,map <unsigned,unsigned char> & Donor_CT,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Acc_GC,map <unsigned,unsigned char> & Acc_AC,map <unsigned,unsigned char> & Acc_AT);
void Match_Unmapped(FILE* Unmapped,Hash & Junctions,COVTYPE* Coverage,map <unsigned,unsigned char> & Acc_AG,map <unsigned,unsigned char> & Acc_AC);
void Match_UnmappedM(FILE* Unmapped,Hash & Junctions,map <unsigned,unsigned char> & Acc_GC,map <unsigned,unsigned char> & Acc_AC,map <unsigned,unsigned char> & Acc_AT);
void Match_UnmappedX(FILE* Unmapped,Hash & Junctions,map <unsigned,unsigned char> & Donor_AT,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Donor_GC);
void Match_UnmappedMX(FILE* Unmapped,Hash & Junctions,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Donor_CT);

char Search_Forwards_Exact(SARANGE & Tag,char* Current_Tag, int Start,int StringLength);
char Extend_Right(SARANGE & Tag,char* Current_Tag,int Start,int StringLength,int Stop_Gap);
char Extend_Left(SARANGE & Tag,char* Current_Tag,int Start,int StringLength,int Stop_Gap);
char Search_Backwards_Exact(SARANGE & Tag,char* Current_Tag,int Start,int StringLength);
char Scan_Half_ExtendR(char* Current_Tag, ORIENTATION & O,unsigned Read_ID);
char Scan_Half_ExtendL(char* Current_Tag, ORIENTATION & O,unsigned Read_ID);
char Scan_Half_ExtendLX(char* Current_Tag, ORIENTATION & O,unsigned Read_ID);
char Guess_Orientation(int Off);
char Load_Tail (char* Current_Tag, SARANGE & Tag);
char Load_Head (char* Current_Tag, ORIENTATION & O,SARANGE & Tag);
char Read_Tag(READ & Head,READ & Tail);
unsigned Location(SARANGE & Tag, char Index);
unsigned Get_File_Size(FILE* File);
unsigned Get_Block_Start(unsigned SAValue,unsigned & M);
gzFile File_OpenZ(const char* File_Name,const char* Mode);
BWT* Load_Indexes(char *BWTINDEX,char *OCCFILE, char *SAFILE);
void UnLoad_Indexes();
FILE* File_Open(const char* File_Name,const char* Mode);
unsigned Splicing_Signal(Island* Current_Island, Island* Next_Island);
void Merge_Island(Island* Prev_Island, Island* Current_Island);
bool Loop_For_Junctions(map <unsigned,unsigned char>::iterator &Start, map<unsigned,unsigned char>::iterator & End, unsigned char Sig, int Pairing, int x, int y, char* Read, Hash & Junctions, int Direction);
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

//{---------------------------- GLOBAL VARIABLES -------------------------------------------------
char Hit_Mask=0;
unsigned Sum=0;
unsigned T=0;
char debug;
char const *Int_To_Junc[10];
string Genome_String;
State* Can_Acceptor;
State* Can_Donor;
MMPool* mmPool;
Offset_Record Genome_Offsets[80];
gzFile Input_File;
gzFile Mate_File;
FILE* Junc_Log;
FILE* Input_FileO;
FILE* F;
READ Head, Tail;
RANGEINDEX Range_Index;
ORIENTATION Orientation[4];//orientation info...
BWT* revfmi,*fwfmi;
SARANGE Cache_SF[MAXTAG+1];
SARANGE Cache_SB[MAXTAG+1];
SARANGE Branch_Ranges[4];
SARANGE Temp_Branch_Ranges[4];
SARANGE* Head_Hits_Pos;//[30];
SARANGE* Head_Hits_Neg;//[30];
SARANGE* Tail_Hits_Neg;//[30];
SARANGE* Tail_Hits_Pos;//[30];
SARANGE* Two_Mismatches_At_EndP;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARANGE* Two_Mismatches_At_EndC;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARANGE* Two_Mismatches_At_End_ForwardP;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARANGE* Two_Mismatches_At_End_ForwardC;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARANGE* Mismatches_BackwardP;//stores possible 2-mismatches
SARANGE* Mismatches_BackwardC;//stores possible 2-mismatches
SARANGE* Mismatches_ForwardP;//stores possible 2-mismatches
SARANGE* Mismatches_ForwardC;//stores possible 2-mismatches
SARANGE* FSSStack;
SARANGE* BMStack;
SARANGE* BMHStack;
PAIR* Pairs;
FILE* Original_File;
FILE* SAM;
FILE* BigBed;
FILE* Junctions_File;
FILE* Map_File;
//--------------------------- Range Index Gloals ---------------
SA* SA_Index;
int* SA_Blocks;
unsigned Hits;
char COMPRESS;
//--------------------------- Range Index Gloals ---------------

char First_Pass=TRUE;
char COUNT_LF_MOTIF=TRUE; 
char MAPMODE=TRUE;
char USEREFGENE=FALSE;
char PROCESSMODE=TRUE;
char NORMAL_TAGS; 
char MINIMUMINDEX=FALSE;	
char PAIRING_TYPE;
char FILETYPE;
char BUILDREF=FALSE;
char MINIMIZE_GAP=FALSE;
char WRITE_SPLITREAD=FALSE;
char DUMP_ALL_JUNC=FALSE;
char Half_Intact;
char Char_To_CodeC[256];
char Char_To_Code[256];
char Packed[MAXTAG+1];
char Random_Array[]="tatacgataggacaatgtcttcgaagcccacgcggtaagccggtcattgcggttgtgcgaacactatcagcctcgctgcatggttaccctgggtggataggacgtttgcccgacattttgacacgcataaaaggtctgtagtgggggtggcacaccataaaccctggggcggctccacgatcgtaaaatcctgcgatctg";
char JunctionFile[]="junctions.bed";
char Wiggles[]="coverage.bed";
char MapFile[]="mapfile";
char Patternfile_Count=0;//number of query files (head /tail) 
char Try[4];
int HITS;
int EXON_JOIN_THRESHOLD=10;//500;//76;

int PAIR_LENGTH_RIGHT,PRRH,PRLH;
int PAIR_LENGTH_LEFT,PLRH,PLLH;
int TAG_COPY_LEN;
int STRINGLENGTH,LH,RH,RHQR,RHQL,STRINGLENGTHO;
int MAX_TAGS_TO_PROCESS;
int MIS_IN_INITMAP=2;
int MIS_SPLICE=2;
int SAGAP_CUTOFF_T,SAGAP_CUTOFF_H;
int RQFACTOR=18;//length of indexed rqindex string..
int MAXHITS=1;
int ARRAY_BOUND,END_BOUND;
int Two_Mismatches_At_End_PointerP; 
int Two_Mismatches_At_End_PointerC; 
int Mismatches_Forward_PointerP;
int Mismatches_Forward_PointerC;
int Mismatches_Backward_PointerP;
int Mismatches_Backward_PointerC;
int Two_Mismatches_At_End_Forward_PointerP=0;
int Two_Mismatches_At_End_Forward_PointerC=0;
int Genome_Count;int COUNT=0;int COUNT_IN_EXT=0;//int COUNT=2;int COUNT_IN_EXT=3;
int MIN_SUPPORT=0; //Minimum number of reads needed to cover a base...
int MAX_HITS_TO_STORE=200;//later resolve this and MAXCOUNT conflict..
int MAXCOUNT=MAX_HITS_TO_STORE; //30000  //Maximum number of pairs to enum...



unsigned File_Size;
unsigned Dummy[]={0,0};
unsigned EXONGAP=1000;
unsigned SOURCELENGTH;
unsigned Branch_Characters[4];
unsigned Temp_BC[4];
unsigned Conversion_Factor;
unsigned CONVERSION_FACTOR;//factor for rqindex..
unsigned Offsets[80];


const char* Code_To_Char="acgt";const char* Code_To_CharCAPS="ACGT";const char* Code_To_CharC="tgca";const char* Code_To_CharCAPSC="TGCA";
char* JUNCTIONFILE=JunctionFile;
char* WIGGLEFILE=Wiggles;
char* MAPFILE=MapFile;
char* PATTERNFILE;
char* PATTERNFILE1;
char* BWTFILE ; 
char* OCCFILE ;
char* SAFILE;
char* REVBWTINDEX;
char* REVOCCFILE;
char* REVSAFILE;
char* LOCATIONFILE;
char* INDFILE;
char* BLKFILE;
char* BINFILE;
char* SORTEDRANGEFILE;
char* RANGEFILE;
char* Source;//buffer for command line processing...
char* Current_Tag;
unsigned char* Original_Text;//encoded genome file...
unsigned Orig_Text;

Out_Record* Buffer;
Out_Record Unique;
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------

//{---------------------------- Command Line  -------------------------------------------------
option Long_Options[]=
{
{"help",0,NULL,'h'},
{"query",1,NULL,'q'},
{0,0,0,0}
};

//}---------------------------- Command Line -------------------------------------------------

int main(int argc, char* argv[])
{
	unsigned Total_Hits=0,Tags_Processed=0,Tag_Count=0;
	time_t Start_Time,End_Time;
	unsigned Number_of_Tags=1000,Average_Length;
	unsigned Progress=0;
	char Pass_Count;
	unsigned Read_ID=0;
	char FULLMATCHONLY=FALSE;//TRUE;

	SARANGE Tag;

	Parse_Command_line(argc,argv);
	if(mkdir("Raw_Out",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) && errno != EEXIST) {printf("Cannot create Temp directory..\n");exit(-1);};
	Genome_Count=Load_Locations(Genome_Offsets,Offsets);
	Init();
	Verbose();
	time(&Start_Time);
	time_t Maptime;
	if(MAPMODE)
	{
		printf("======================]\r[");//progress bar....
		if (NORMAL_TAGS) Pass_Count=1; else Pass_Count=2;//Number of passes to make..
		while (Read_Tag(Head,Tail))
		{
			char Head_Tail=0;
			if (MAX_TAGS_TO_PROCESS && Tag_Count >= MAX_TAGS_TO_PROCESS) break; 
			while (Head_Tail<Pass_Count)//loop head, then tail
			{
				if(Head_Tail) Head=Tail;
				Head_Tail++;
				Tag_Count++;
				Progress++;
				if (Progress==Number_of_Tags) 
				{
					if (MAX_TAGS_TO_PROCESS)
					{
						Number_of_Tags=(MAX_TAGS_TO_PROCESS)/20;
						Progress=0;
						Show_Progress(Tag_Count*100/MAX_TAGS_TO_PROCESS);
					}
					else
					{
						off64_t Current_Pos=ftello64(Input_FileO);
						Average_Length=Current_Pos/Tag_Count+1;//+1 avoids divide by zero..
						Number_of_Tags=(File_Size/Average_Length)/20;
						Progress=0;
						Show_Progress(Current_Pos*100/File_Size);
					}
				}

				Read_ID++;
				if(!Scan(Total_Hits,Tags_Processed,0,Read_ID))
				{
					int olds=STRINGLENGTH;
					int oldi=MIS_IN_INITMAP;
					int oldc=COUNT_IN_EXT;
					MIS_IN_INITMAP=0;
					STRINGLENGTH=50;
					//COUNT_IN_EXT=1;
					Scan(Total_Hits,Tags_Processed,0,Read_ID);
					//Scan(Total_Hits,Tags_Processed,0,Read_ID);
					{
						Read_ID++;
						int oldso=STRINGLENGTHO;
						int o=STRINGLENGTHO-50;
						STRINGLENGTHO=STRINGLENGTH=50;
						Scan(Total_Hits,Tags_Processed,o,Read_ID);
						//Scan(Total_Hits,Tags_Processed,o,Read_ID);
						{
							Read_ID++;
							STRINGLENGTHO=oldso;
							Scan(Total_Hits,Tags_Processed,18,Read_ID);//19
						}
						STRINGLENGTHO=oldso;
//left split..
					}
					STRINGLENGTH=olds;
					MIS_IN_INITMAP=oldi;
					COUNT_IN_EXT=oldc;
				}
			}
		}
		UnLoad_Indexes();
		printf("\r[++++++++100%%+++++++++]\n");//progress bar....
		time(&End_Time);
		Maptime=difftime(End_Time,Start_Time);
	}

	if (PROCESSMODE)
	{
		for(int i=0;i<Genome_Count;i++)
		{
			Genome_Offsets[i].Offset=Genome_Offsets[i+1].Offset;
			printf("Processing Chromosome %s...\n",Genome_Offsets[i].Genome);
			Orig_Text= Offsets[i];
			Coverage(Genome_Offsets[i]);
		}
	}
	printf("Average %d\n",Sum);
	printf("%u / %u Tags/Hits",Tag_Count,Tags_Processed+Total_Hits);
	time(&End_Time);printf("\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	printf("\n Time Taken  - %.0lf Seconds ..\n ",Maptime);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Guess_Orientation
 *  Description:  Try to guess the strand of the read...
 *		  Scanning is started from Read[Off,...,Off+STRINGLENGTH]
 *  		  return 0 for head, 1 for tail
 *  		  if TAGRQF.Start=0, read mapped ...
 * =====================================================================================
 */
char Guess_Orientation(int Off)
{
	char Half_OK_Read=FALSE;
	char Half_OK_Compl=TRUE;

	SARANGE TagF,TagFR,TagFH;//Tag,Tag for RQFactor, Tag for LH.
	SARANGE TagFC,TagFCR,TagFCH;
	SARANGE TagB,TagBR,TagBH;//Tag,Tag for RQFactor, Tag for LH.
	SARANGE TagBC,TagBCR,TagBCH;

	Mismatches_Forward_PointerP=0;//init pointers..
	Mismatches_Forward_PointerC=0;
	Mismatches_Backward_PointerP=0;
	Mismatches_Backward_PointerC=0;
	Two_Mismatches_At_End_PointerP=0;
	Two_Mismatches_At_End_PointerC=0;
	Two_Mismatches_At_End_Forward_PointerP=0;
	Two_Mismatches_At_End_Forward_PointerC=0;

	HITS=0;Orientation[0].Strand=0;
	TagF.Start=1;TagF.End=SOURCELENGTH;TagF.Skip=0;TagF.Mismatches=0;TagF.Mismatch_Char=0;
	TagFC.Level=TagB.Level=TagBC.Level=0;
	TagFC.Start=TagB.Start=TagBC.Start=0;
	Try[LEFT]=0;Try[RIGHT]=0;
	Try[LEFTC]=0;Try[RIGHTC]=0;

	Current_Tag=Head.Tag+Off;
	Cache_SF[RQFACTOR].Start=0;
	Cache_SF[LH].Start=0;
	char zfound=0,zcfound=0;Hit_Mask=0;
	if((zfound=Search_Forwards_Exact(TagF,Head.Tag+Off,1,STRINGLENGTH)) )//search + for exact..
	{
		if(TagF.Start==TagF.End || TagF.Skip)//Unique hit..
		{
			Print_LocationX(TagF,FW);
			HITS++;
			//return TRUE;//exact match...
		}
		else
		{
			HITS+=TagF.End-TagF.Start+1;
		}
		//else
		{
			//return TRUE;//exact match...
			TagFR=Cache_SF[RQFACTOR];//save important sa ranges..
			TagFH=Cache_SF[LH];
			Try[LEFT]=1;//check feasibility of extentions..

			TagB.Start=1;TagB.End=SOURCELENGTH;TagB.Skip=0;TagB.Mismatches=0;TagB.Mismatch_Char=0;
			Search_Backwards_Exact(TagB,Head.Tag+Off,STRINGLENGTH,STRINGLENGTH);

			Try[RIGHT]=1;
			TagBH=Cache_SB[RH];
			TagBR=Cache_SB[RQFACTOR];
			Hit_Mask=TRUE;
		}
	}
	//else  //search - for exact...
	{

		TagFR=Cache_SF[RQFACTOR];//save important sa ranges..
		TagFH=Cache_SF[LH];
		Cache_SF[RQFACTOR].Start=0;
		Cache_SF[LH].Start=0;
		Try[LEFT]=(TagF.Level >= LH);//check feasibility of extentions..

		TagFC.Start=1;TagFC.End=SOURCELENGTH;TagFC.Skip=0;TagFC.Mismatches=0;TagFC.Mismatch_Char=0;Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off);
		if((zcfound=Search_Forwards_Exact(TagFC,Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off),1,STRINGLENGTH)))//search complement for exact...
		{
			if(TagFC.Start==TagFC.End || TagFC.Skip)//Unique hit..
			{
				Print_LocationX(TagFC,FW);
				HITS++;
				//return TRUE;//exact match..
			}
			else
			{
				HITS+=TagFC.End-TagFC.Start+1;
			}
			//else
			{
				TagFCR=Cache_SF[RQFACTOR];//save important sa ranges..
				TagFCH=Cache_SF[LH];
				Try[LEFTC]=1;//check feasibility of extentions..

				TagBC.Start=1;TagBC.End=SOURCELENGTH;TagBC.Skip=0;TagBC.Mismatches=0;TagBC.Mismatch_Char=0;
				Search_Backwards_Exact(TagBC,Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off),STRINGLENGTH,RH+1); 
				Try[RIGHTC]=1;
				TagBCH=Cache_SB[RH];
				TagBCR=Cache_SB[RQFACTOR];

				Hit_Mask=TRUE;
			}
		}
		//if(!Hit_Mask)
		{
			TagFCR=Cache_SF[RQFACTOR];//save important sa ranges..
			TagFCH=Cache_SF[LH];
			Try[LEFT]=(TagF.Level >= LH);//check feasibility of extentions..
			Try[LEFTC]=(TagFC.Level >=LH);

			if ( MIS_IN_INITMAP && (Try[LEFT] || Try[LEFTC]) && !(zcfound ||zfound) )//left half intact..
			{
				if (TagF.Level > TagFC.Level)//guess orientation..
				{
					SARANGE Tag=TagFH;Current_Tag=Head.Tag+Off;Tag.Level=1;
					if(Search_Forwards(Tag,1,LH+1,RH,Two_Mismatches_At_End_Forward_PointerP,Two_Mismatches_At_End_ForwardP,Mismatches_Forward_PointerP,Mismatches_ForwardP)) HITS++;//return TRUE;

					if (Try[LEFTC])//complement is also long..
					{
						SARANGE Tag=TagFCH;Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off);Tag.Level=1;
						if(Search_Forwards(Tag,1,LH+1,RH,Two_Mismatches_At_End_Forward_PointerC,Two_Mismatches_At_End_ForwardC,Mismatches_Forward_PointerC,Mismatches_ForwardC)) HITS++;//return TRUE;
					}
					else {Half_Intact=FALSE;Orientation[1].Strand=0;}//there are breaks in both halves..
				}
				else
				{
					SARANGE Tag=TagFCH;Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off);Tag.Level=1;
					if(Search_Forwards(Tag,1,LH+1,RH,Two_Mismatches_At_End_Forward_PointerC,Two_Mismatches_At_End_ForwardC,Mismatches_Forward_PointerC,Mismatches_ForwardC)) HITS++;//return TRUE;

					if (Try[LEFT])//complement is also long..
					{
						SARANGE Tag=TagFH;Current_Tag=Head.Tag+Off;Tag.Level=1;
						if(Search_Forwards(Tag,1,LH+1,RH,Two_Mismatches_At_End_Forward_PointerP,Two_Mismatches_At_End_ForwardP,Mismatches_Forward_PointerP,Mismatches_ForwardP)) HITS++;// return TRUE;
					}
					else {Half_Intact=FALSE;Orientation[1].Strand=0;}//There are breaks in both halves..
				}

			}
		}
	}
	//else//break in left half..
	if(1)// && !(zfound || zcfound)) 
	{
		TagB.Start=1;TagB.End=SOURCELENGTH;TagB.Skip=0;TagB.Mismatches=0;TagB.Mismatch_Char=0;
		Search_Backwards_Exact(TagB,Head.Tag+Off,STRINGLENGTH,STRINGLENGTH);

		if (Try[RIGHT]=(TagB.Level >= RH))
		{
			TagBH=Cache_SB[RH];
			TagBR=Cache_SB[RQFACTOR];

			//do one mismatch...
			if (MIS_IN_INITMAP)
			{
				SARANGE Tag=TagBH;Current_Tag=Head.Tag+Off;Tag.Level=1;
				if (Search_Backwards(Tag,1,LH,LH,Two_Mismatches_At_End_PointerP,Two_Mismatches_At_EndP,Mismatches_Backward_PointerP,Mismatches_BackwardP)) HITS++;//return TRUE;
			}
		}

		TagBC.Start=1;TagBC.End=SOURCELENGTH;TagBC.Skip=0;TagBC.Mismatches=0;TagBC.Mismatch_Char=0;
		Search_Backwards_Exact(TagBC,Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off),STRINGLENGTH,STRINGLENGTH); 
		if (Try[RIGHTC]=(TagBC.Level >= RH))
		{
			TagBCH=Cache_SB[RH];
			TagBCR=Cache_SB[RQFACTOR];

			//do one mismatch...
			if (MIS_IN_INITMAP)
			{
				SARANGE Tag=TagBCH;Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off);Tag.Level=1;
				if (Search_Backwards(Tag,1,LH,LH,Two_Mismatches_At_End_PointerC,Two_Mismatches_At_EndC,Mismatches_Backward_PointerC,Mismatches_BackwardC)) HITS++;//return TRUE;
			}
		}

	}
//---------------------------------------------------------------------------------------------------------------------------------------------------			
        //do two mismatch 2|0..
	//int MaxF=(TagF.Level>TagFC.Level) ? TagF.Level : TagFC.Level;
	//int MaxB=(TagB.Level>TagBC.Level) ? TagB.Level : TagBC.Level;
	int MaxF=(TagF.Level>TagB.Level) ? TagF.Level : TagB.Level;
	int MaxB=(TagBC.Level>TagFC.Level) ? TagBC.Level : TagFC.Level;
	if(MaxF>MaxB)//guess orientation...
	{
		Orientation[0].Strand='+';
		Orientation[0].String=Head.Tag+Off;
		Orientation[0].TagF=TagF;
		Orientation[0].TagFR=TagFR;
		Orientation[0].TagFH=TagFH;
		Orientation[0].TagB=TagB;
		Orientation[0].TagBR=TagBR;
		Orientation[0].TagBH=TagBH;

		Orientation[1].Strand='-';
		Orientation[1].String=Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off);
		Orientation[1].TagF=TagFC;
		Orientation[1].TagFR=TagFCR;
		Orientation[1].TagFH=TagFCH;
		Orientation[1].TagB=TagBC;
		Orientation[1].TagBR=TagBCR;
		Orientation[1].TagBH=TagBCH;
	}
	else 
	{
		Orientation[1].Strand='+';
		Orientation[1].String=Head.Tag+Off;
		Orientation[1].TagF=TagF;
		Orientation[1].TagFR=TagFR;
		Orientation[1].TagFH=TagFH;
		Orientation[1].TagB=TagB;
		Orientation[1].TagBR=TagBR;
		Orientation[1].TagBH=TagBH;

		Orientation[0].Strand='-';
		Orientation[0].String=Head.Complement+(STRINGLENGTHO-STRINGLENGTH-Off);
		Orientation[0].TagF=TagFC;
		Orientation[0].TagFR=TagFCR;
		Orientation[0].TagFH=TagFCH;
		Orientation[0].TagB=TagBC;
		Orientation[0].TagBR=TagBCR;
		Orientation[0].TagBH=TagBCH;
	}
	if (HITS) {return TRUE;}
	if (MIS_IN_INITMAP < 2 ) return FALSE;
	if(!(zcfound || zfound)) return Two_Mismatch_Search();
//---------------------------------------------------------------------------------------------------------------------------------------------------			
	//Should have found all exact, 1 mismatches and 2|0 mismatches...
}
//{--------------------------------------  Coverage -------------------------------------------------------------------
void Coverage(Offset_Record Current_Genome)
{
	Hash Junctions;
	Hash Unmapped_Junctions;
	OP JPair;
	OPX* Junc_List=new OPX[MAX_JUNC_LIST+1];
	int Junc_List_Ptr=MAX_JUNC_LIST;
	unsigned Last_ID=0;

	COVTYPE* Coverage;

	int Islands=0;
	int Max_Islands= -1;//number of nodes in the linked list


	if(!(Coverage=(COVTYPE*)calloc(Current_Genome.Offset,sizeof(COVTYPE)))) {printf("Coverage():Malloc error!..\n");exit(0);}
	//if(!(Confident_Coverage=(COVTYPE*)calloc(Current_Genome.Offset,sizeof(COVTYPE)))) {printf("Coverage():Malloc error!..\n");exit(0);}
	char Pass=0,Strand=PLUS;
	char Parsed=FALSE;

	while(Pass<2)
	{
		FILE* Handle=(Strand == PLUS ) ? Current_Genome.Out_File : Current_Genome.Out_FileM;
		rewind(Handle);
		char E_O_F=FALSE;
		while (!E_O_F)
		{
			int Chars_In_Buffer=fread(Buffer,sizeof(Out_Record),BUFSIZE,Handle);
			if (Chars_In_Buffer<BUFSIZE) E_O_F=TRUE;
			int Current_Junc_Pos=0;
			while(Current_Junc_Pos < Chars_In_Buffer)
			{
				STRINGLENGTH=Buffer[Current_Junc_Pos].Strlen;
				unsigned Tail_Start=Buffer[Current_Junc_Pos].TLocation;//tail portion..
				int Tail_Length=Buffer[Current_Junc_Pos].TLength;
				int Tail_Stop=Tail_Start+Tail_Length;
				unsigned Head_Start=Buffer[Current_Junc_Pos].HLocation;//Head portion..
				int Head_Length=Buffer[Current_Junc_Pos].HLength;
				int Head_Stop=Head_Start+Head_Length;
				//int Junc_Count=Buffer[Current_Junc_Pos].Junc_Count;
				char Type=Buffer[Current_Junc_Pos].Type;
				unsigned ID=Buffer[Current_Junc_Pos].ID;
				Current_Junc_Pos++;
				Parsed=TRUE;
//Parse junctions..
				if (Type != FULLREAD && Last_ID != ID) //New read..
				{
					if (Junc_List_Ptr < MAX_JUNC_LIST)//not too many juncs or init round..
					{
						Find_Split_Junctions(Junc_List,Junc_List_Ptr,Junctions,Last_ID,Current_Genome.Genome);
					}
					Junc_List_Ptr=0;
					Last_ID=ID;
				}
//Build coverage ...........
				assert(Tail_Stop<Current_Genome.Offset && Head_Stop<Current_Genome.Offset && Head_Start<Head_Stop && ((Tail_Start==0 && Tail_Stop==0)||(Tail_Start<Tail_Stop)));
				for (int i=Tail_Start;i<Tail_Stop;i++) if(Coverage[i] < COVTYPEMAX) Coverage[i]++;
				for (int i=Head_Start;i<Head_Stop;i++) if(Coverage[i] < COVTYPEMAX) Coverage[i]++;

//Collect possible junctions ...........
				if (Type !=FULLREAD && Junc_List_Ptr<MAX_JUNC_LIST)
				{
					if(Head_Length && Tail_Length && Head_Length + Tail_Length >= STRINGLENGTH)//both parts have been mapped...
					{
						if (Head_Start+STRINGLENGTH > Tail_Start) //What is the purpose of this step?
						{
							Junc_List[Junc_List_Ptr].x=Head_Start;
							Junc_List[Junc_List_Ptr].y=Tail_Start;
							Junc_List[Junc_List_Ptr++].Motif=0;
						} 
						else Juncs_In_Read(Junc_List,Junc_List_Ptr,Head_Start,Head_Length,Tail_Start,Tail_Length);
					}
					else//Not functional yet.. need to save to same file..
					{
						Junc_List[Junc_List_Ptr].x=Head_Start;
						Junc_List[Junc_List_Ptr].y=Tail_Start;
						Junc_List[Junc_List_Ptr++].Motif=0;
					}
				}
			}
		}
		Pass++;Strand=MINUS;
	}
	if (!Parsed) {free(Coverage);delete [] Junc_List;return;}//no coverage...
	//scan islands and find transcriptome...
	char In_Island=FALSE;int ScanPos=0;
	Island* Island_List=NULL;Island* Current_Island=NULL;
	unsigned Last;

	for (;ScanPos<Current_Genome.Offset;ScanPos++)
	{
		if(Coverage[ScanPos] < MIN_SUPPORT) Coverage[ScanPos]=0;//Filter out low coverage reads... MIN_SUPPORT =0 if no filtering...
		if (Coverage[ScanPos])
		{
			if (!In_Island)//starting an island..
			{
				//printf("[%u",ScanPos);
				if (Current_Island)
				{
					Current_Island->Next=new Island;
					Current_Island=Current_Island->Next;
					
				}
				else {Current_Island=new Island;Island_List=Current_Island;}
				Current_Island->Next= NULL;
				Current_Island->Status=0;
				Current_Island->Start=ScanPos;
				In_Island=TRUE;
			}
		}
		else if (In_Island)
		{
			Current_Island->End=ScanPos;
			Islands++;
			//printf(":%u]\n",ScanPos);
			In_Island=FALSE;
		}
	}
	if (In_Island) {Current_Island->End=ScanPos;}//printf(":%u]\n",ScanPos);}
	Max_Islands = (Max_Islands > Islands) ? Islands : Max_Islands;


    //======================Merging Islands=====================
    if(Island_List && Island_List->Next)
    {
        Island* My_Prev_Island = Island_List;
        Island* My_Current_Island = Island_List->Next;
        unsigned My_Prev_Junc_Signal = Splicing_Signal(My_Prev_Island, My_Current_Island);
        unsigned My_Current_Junc_Signal = Splicing_Signal(My_Current_Island, My_Current_Island->Next);
        unsigned My_Next_Junc_Signal;
        while(My_Current_Island->Next) {
            Island* My_Next_Island = My_Current_Island->Next;
            if(My_Current_Junc_Signal > 0 && My_Prev_Junc_Signal == 0
                    && My_Current_Island->Start - My_Prev_Island->End < 10)
            {
                Merge_Island(My_Prev_Island, My_Current_Island);
            }
            My_Next_Junc_Signal = Splicing_Signal(My_Next_Island,My_Next_Island->Next);
            if(My_Current_Junc_Signal > 0 && My_Next_Island->Next && My_Next_Junc_Signal == 0)
            {
                Merge_Island(My_Next_Island,My_Next_Island->Next);
            }
            My_Current_Island = My_Current_Island->Next;
        }
    }
	printf("=======================================================\n");

	Hash RefGene;
	JUNC_HASH Donor_AT,Donor_GT,Donor_GC,Donor_CT;
	JUNC_HASH Acc_GC,Acc_AC,Acc_AT,Acc_AG;
	if (USEREFGENE) 
	{
		FILE* Ref_Handle=Current_Genome.Ref_File; 
		Load_RefGene(Ref_Handle,RefGene,Donor_AT,Donor_GT,Donor_GC,Donor_CT,Acc_AC,Acc_AG,Acc_AT,Acc_GC);
	}
//---------------------------------------------------------------------------------------------------------------------------

	Seek_Island_Junctions(Island_List,Donor_AT,Donor_GT,Donor_GC,Acc_AG,Acc_AC);
	Seek_Island_JunctionsM(Island_List,Donor_CT,Donor_GT,Acc_GC,Acc_AC,Acc_AT);

//------------------------------- Map unmapped Left half long --------------------------------
	FILE* Unmapped=Current_Genome.Unmapped; 
	Match_Unmapped(Unmapped,Junctions,Coverage,Acc_AG,Acc_AC);
	Match_UnmappedM(Unmapped,Junctions,Acc_GC,Acc_AC,Acc_AT);

    printf("After Match_Unmapped Unmapped\n");

	Unmapped= Current_Genome.UnmappedM;
	Match_Unmapped(Unmapped,Junctions,Coverage,Acc_AG,Acc_AC);
	Match_UnmappedM(Unmapped,Junctions,Acc_GC,Acc_AC,Acc_AT);
	
//------------------------------- Map unmapped Right half long --------------------------------
	FILE* UnmappedX=Current_Genome.UnmappedX; 
	Match_UnmappedX(UnmappedX,Junctions,Donor_AT,Donor_GT,Donor_GC);
	Match_UnmappedMX(UnmappedX,Junctions,Donor_GT,Donor_CT);

	UnmappedX=Current_Genome.UnmappedXM; 
	Match_UnmappedX(UnmappedX,Junctions,Donor_AT,Donor_GT,Donor_GC);
	Match_UnmappedMX(UnmappedX,Junctions,Donor_GT,Donor_CT);

//--------------------------------------------- Enter Exons to Hash ----------------------------
// These are the putative exons from coverage...
	HASH_PAIR ExonH;
	HASH_PAIR ExonT;
	HASH_PAIR::iterator ExonH_I;
	HASH_PAIR::iterator ExonT_I;
	Current_Island=Island_List;
	while(Current_Island)
	{
		if (Current_Island->Status != DELETED)
		{
			ExonH[Current_Island->Start]=Current_Island->End;//ExonH <start exon, end exon>
			ExonT[Current_Island->End]=Current_Island->Start;//ExonT <end exon,start exon>
		}
		Current_Island=Current_Island->Next;
	}
	HASH_PAIR Acceptor;
	HASH_PAIR Donor;
	JUNC_HASH::iterator Junc_I;
	for (Junc_I=Acc_AG.begin();Junc_I != Acc_AG.end();Junc_I++) Acceptor[Junc_I->first]=S_ag;
	for (Junc_I=Acc_AT.begin();Junc_I != Acc_AT.end();Junc_I++) Acceptor[Junc_I->first]=S_at;
	for (Junc_I=Acc_GC.begin();Junc_I != Acc_GC.end();Junc_I++) Acceptor[Junc_I->first]=S_gc;
	for (Junc_I=Acc_AC.begin();Junc_I != Acc_AC.end();Junc_I++) Acceptor[Junc_I->first]=S_ac;
    //for (Junc_I=Acc_GT.begin();Junc_I != Acc_GT.end();Junc_I++) Acceptor[Junc_I->first]=S_gt;

	for (Junc_I=Donor_CT.begin();Junc_I != Donor_CT.end();Junc_I++) Acceptor[Junc_I->first]=S_ct;
	for (Junc_I=Donor_AT.begin();Junc_I != Donor_AT.end();Junc_I++) Acceptor[Junc_I->first]=S_at;
	//for (Junc_I=Donor_AC.begin();Junc_I != Donor_AC.end();Junc_I++) Acceptor[Junc_I->first]=S_ac;
	for (Junc_I=Donor_GT.begin();Junc_I != Donor_GT.end();Junc_I++) Acceptor[Junc_I->first]=S_gt;
	for (Junc_I=Donor_GC.begin();Junc_I != Donor_GC.end();Junc_I++) Acceptor[Junc_I->first]=S_gc;

	//for (Exons_I=Exons.begin();Exons_I != Exons.end();Exons_I++) printf ("{%d:%d}\n",Exons_I->first,Exons_I->second);
//--------------------------------------------- Enter Exons to Hash --------------------------------------------------------
	Print_BED(Junctions,Island_List,Coverage,Current_Genome.Genome,Strand);

	Current_Island=Island_List;//free memory..
	while (Current_Island)
	{
		Island* Next=Current_Island->Next;
		delete Current_Island;
		Current_Island=Next;
	}
	free (Coverage);
	delete []Junc_List;

}

void Find_Split_Junctions(OPX* Junc_List,int Junc_List_Ptr,Hash & Junctions,unsigned Read_ID,char* Genome)
{
#define MAX_JUNC 50
#define dist(x,y) (((x)>(y)) ? (x)-(y): (y)-(x))


	OP HPair,LPair;
	int Hi_Count=0,LF_Motif_Count=0,Motif,MotifH,MotifL,J_Ptr=0;
	unsigned x,y;

	//if (Junc_List_Ptr >1) return;
	for (int i=Junc_List_Ptr-1;i>=0;i--)
	{

		Motif=Junc_List[i].Motif;
		x=Junc_List[i].x;
		y=Junc_List[i].y;
		if (Motif == S_gt_ag ||Motif == S_ct_ac)
		{
			//if (!(Hi_Count && HPair.x==x && HPair.y==y))//not a duplicate junction...
			{
				bool Reject=false;
				for (int j=0;j<=Junc_List_Ptr-1;j++)
				{
					int MotifT=Junc_List[j].Motif;
					if ((i !=j) && (MotifT == S_gt_ag ||MotifT == S_ct_ac||MotifT==0))
					{
						if (Junc_List[j].x == x)
						{
							if (Junc_List[j].y==y){Reject=true;Junc_List[(i<j? i:j)].Motif=0;}//Duplicate junction..
							else if (dist(Junc_List[j].y,y)<400000){Reject=true;Junc_List[i].Motif=Junc_List[j].Motif=0;}//Duplicate junction..

						}
						else
						{
							//assert (Junc_Found[j].y == y);
							if (dist(Junc_List[j].x,x)<400000){Reject=true;Junc_List[i].Motif=Junc_List[j].Motif=0;}//Duplicate junction..
						}
					}
				}
			}
		}
		/*else
		{
			if (!(LF_Motif_Count && LPair.x==x && LPair.y==y))//not a duplicate junction...
			{
				LPair.x=x;
				LPair.y=y;
				MotifL=Motif;
				LF_Motif_Count++;
			}
		}*/
	}

	int J_Major_Splice=0,J_Minor_Splice=0;//Are there multuple juncs?
	for (int i=0;i<Junc_List_Ptr;i++)
	{
		int Motif=Junc_List[i].Motif;
		if (Motif)
		{
			if (Motif == S_gt_ag ||Motif == S_ct_ac) 
			{
				J_Major_Splice++; 
				HPair.x=Junc_List[i].x;
				HPair.y=Junc_List[i].y;
				MotifH=Junc_List[i].Motif;
			}
			else 
			{
				J_Minor_Splice++; 
				LPair.x=Junc_List[i].x;
				LPair.y=Junc_List[i].y;
				MotifL=Junc_List[i].Motif;
			}
		}
	}
	//if (J_Major_Splice+J_Minor_Splice==1) Junctions.Insert(HPair,MotifH,10)//strict..
	if (J_Major_Splice==1) {Junctions.Insert(HPair,MotifH);fprintf(Junc_Log,"%s\t%u\t%u\t%u\n",Genome,HPair.x-2,HPair.y+2+1,Read_ID);}
	else if (J_Major_Splice+J_Minor_Splice==1) {Junctions.Insert(LPair,MotifL);fprintf(Junc_Log,"%s\t%u\t%u\t%u\n",Genome,LPair.x-2,LPair.y+2+1,Read_ID);}

	

}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Mark_Exons
 *  Description:  Convert a list of exons in to ordered set of exons...
 *  		  Use only in index building step...
 * =====================================================================================
 */
void Mark_Exons()
{

	Xon Exon;
	Buffer=(Out_Record*)malloc(BUFSIZE*sizeof(Out_Record));

	Genome_Count=Load_Locations(Genome_Offsets,Offsets);
	printf("Processing refgene....\n");
	while (--Genome_Count >=0)
	{
		for (char PM=0;PM<2;PM++)
		{
			Offset_Record Current_Genome=Genome_Offsets[Genome_Count];
			if (PM) {sprintf((char*)Buffer+sprintf((char*)Buffer,"%s",Current_Genome.Genome),".-.bin");}
			else {sprintf((char*)Buffer+sprintf((char*)Buffer,"%s",Current_Genome.Genome),".+.bin");}
			printf("Converting %s...",(char*)Buffer);
			FILE* Bin_Plus=File_Open((char*)Buffer,"r+b");
			unsigned GSize=Genome_Offsets[Genome_Count+1].Offset;
			unsigned char* Coverage=(unsigned char*)calloc(GSize,sizeof(unsigned char));
			while(fread(&Exon,sizeof(Exon),1,Bin_Plus))//generate coverage..
			{
				for (int i=Exon.Start;i<Exon.End;i++) 
				{
					Coverage[i]++;
				}
			}

			char In_Island=FALSE;
			rewind(Bin_Plus);
			for (unsigned i=0;i<GSize;i++)
			{
				if (Coverage[i])
				{
					if (!In_Island)
					{
						fwrite(&i,sizeof(unsigned),1,Bin_Plus);
						In_Island=TRUE;
					}
				}
				else if (In_Island)
				{
					fwrite(&i,sizeof(unsigned),1,Bin_Plus);
					In_Island=FALSE;
				}
			}
			free (Coverage);
			fclose(Bin_Plus);printf("Done...\n");
		}
	}
}

#define OVERHANG	10

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Scan_Minus_Motif( x, y)
 *  Description:  Check for - strand motif in junction x,y
 * =====================================================================================
 */
int Scan_Minus_Motif(unsigned x,unsigned y)
{
	unsigned char DonorA= (unsigned char)(*(Original_Text+(Orig_Text+x+1)/4) << (((Orig_Text+x+1) % 4) * 2))>>6;
	unsigned char DonorB= (unsigned char)(*(Original_Text+(Orig_Text+x+2)/4) << (((Orig_Text+x+2) % 4) * 2))>>6;

	unsigned char AccA= (unsigned char)(*(Original_Text+(Orig_Text+y-2)/4) << (((Orig_Text+y-2) % 4) * 2))>>6;
	unsigned char AccB= (unsigned char)(*(Original_Text+(Orig_Text+y-1)/4) << (((Orig_Text+y-1) % 4) * 2))>>6;

	int Pairing =0;
	//------------------------------ Enumerate Donor sites ------------------------------------------------------
	// [GT-AG]-> ct - ac
	// [GC-AG]-> ct - gc
	// [AT-AC]-> gt - at
	if (DonorB == N_t)
	{
		if(DonorA == N_c)//CT
			Pairing = S_ct;
		else if (DonorA == N_g)//AT
			Pairing = S_gt;
	}
	//------------------------------ Enumerate Acceptor sites ------------------------------------------------------
	//------------------------------ Enumerate Donor sites ------------------------------------------------------
	if (Pairing == S_ct)
	{
		Pairing=0;
		if (AccB == N_c)
		{
			if(AccA == N_a)//AC [CA-TC]
				Pairing= S_ct_ac;
			else if(AccA == N_g)//CG [CG-TC]
				Pairing= S_ct_gc;
		} 
	}
	else if (Pairing == S_gt) 
	{
		if (AccA == N_a && AccB == N_t)//AT [TA-TG]
			Pairing = S_gt_at; 
		else
			Pairing = 0;
	}

	if (DUMP_ALL_JUNC && !Pairing) Pairing = S_dummy;//dont filter juncs...
	return Pairing;

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Scan_Plus_Motif( x, y)
 *  Description:  Check for + strand motif in junction x,y
 * =====================================================================================
 */
int Scan_Plus_Motif(unsigned x,unsigned y)
{
	unsigned char DonorA= (unsigned char)(*(Original_Text+(Orig_Text+x+1)/4) << (((Orig_Text+x+1) % 4) * 2))>>6;
	unsigned char DonorB= (unsigned char)(*(Original_Text+(Orig_Text+x+2)/4) << (((Orig_Text+x+2) % 4) * 2))>>6;

	unsigned char AccA= (unsigned char)(*(Original_Text+(Orig_Text+y-2)/4) << (((Orig_Text+y-2) % 4) * 2))>>6;
	unsigned char AccB= (unsigned char)(*(Original_Text+(Orig_Text+y-1)/4) << (((Orig_Text+y-1) % 4) * 2))>>6;


	int Pairing =0;
	//------------------------------ Enumerate Acceptor sites ------------------------------------------------------
	if (AccA == N_a)
	{
		if(AccB == N_g)//AG -high chance case [GT-AG] or [GC - AG]
			Pairing = S_ag;
		else if (AccB == N_c)//AC - [AT-AC]
			Pairing = S_ac;
	}
	//------------------------------ Enumerate Acceptor sites ------------------------------------------------------
	//------------------------------ Enumerate Donor sites ------------------------------------------------------
	if (Pairing == S_ag)
	{
		Pairing=0;
		if (DonorA == N_g)
		{
			if(DonorB == N_t)//GT [GT-AG]
				Pairing= S_gt_ag;
			else if(DonorB == N_c)//GC [GC-AG]
				Pairing= S_gc_ag;
		} 
	}
	else if (Pairing == S_ac) 
	{
		if (DonorA == N_a && DonorB == N_t)//AT [AT-AC]
			Pairing = S_at_ac; 
		else
			Pairing = 0;
	}
	return Pairing;
}

#define OVERHANG	10

#define Bit2Char(y) (unsigned char)(*(Original_Text+(Orig_Text+(y))/4) << (((Orig_Text+(y)) % 4) * 2))>>6

void Merge_Island(Island* Prev_Island, Island* Current_Island) {
    Prev_Island->End = Current_Island->End;
    Current_Island->Status = DELETED;
}

unsigned Splicing_Signal(Island* Current_Island, Island* Next_Island) {
    int Pairing [3]= {0,0,0}; //AT-AC,GT-AG,GC-AG
    int Acc [2]= {0,0}; //AC,AG
    int Donor [3]= {0,0,0}; //AT,GT,GC
    int AccM [3] = {0,0,0}; //AT,AC,GC
    int DonorM [2] = {0,0}; //GT,CT
    int PairingM [3] = {0,0,0}; //GT-AT,CT-AC,CT-GC

    if(!Current_Island || !Next_Island)
        return -2;
    else if(Current_Island->Status == DELETED || Next_Island->Status == DELETED)
        return -1;
    
    if(Current_Island && Next_Island && Current_Island->Status != DELETED && Next_Island->Status != DELETED) 
    {
        for (unsigned y=Next_Island->Start+EXTENDISLAND;y<Next_Island->Start-EXTENDISLAND;y--)
        {
            unsigned char AccA= Bit2Char(y-2);
            unsigned char AccB= Bit2Char(y-1);

            if (AccA == N_a) 
            {
                if(AccB == N_g)
                {  
					Acc[1] += 1;
                }
                else if (AccB == N_c)
                {
                    Acc[0] += 1;
                    AccM[1] += 1;
                }
                else if (AccB == N_t)
                {
                    AccM[0] += 1;
                }
            }
            if (AccA == N_g && AccB == N_c)
            {
                AccM[2] += 1;
            }
        }

        for (unsigned x=Current_Island->End-EXTENDISLAND;x>Current_Island->End+EXTENDISLAND;x++)
        {
            unsigned char DonorA= Bit2Char(x+1);
            unsigned char DonorB= Bit2Char(x+2);
            if (DonorA == N_g)
            {
                if(DonorB == N_t)
                {
                    Donor[1] += 1;
                    DonorM[0] += 1;
                }
                else if(DonorB == N_c)
                {
                    Donor[2] += 1;
                }
            }
            if(DonorA == N_a && DonorB == N_t)
            {
                Donor[0] += 1;
            }
            if(DonorA == N_c && DonorB == N_t)
            {
                DonorM[1] += 1;
            }
        }
        if (Donor[0] && Acc[0]) Pairing[0] = S_at_ac;
        if (Donor[1] && Acc[1]) Pairing[1] = S_gt_ag;
        if (Donor[2] && Acc[1]) Pairing[2] = S_gc_ag;
        if (DonorM[0] && AccM[0]) PairingM[0] = S_gt_at;
        if (DonorM[1] && AccM[1]) PairingM[1] = S_ct_ac;
        if (DonorM[1] && AccM[2]) PairingM[2] = S_ct_gc;
    }
    if (Pairing[0] + Pairing[1] + Pairing[2] ||
            PairingM[0] + PairingM[1] + PairingM[2])
        return 2;
    if ((Acc[0] + Acc[1] && Donor[0]+Donor[1]+Donor[2]) ||
            (AccM[0] + AccM[1]+AccM[2] && DonorM[0]+DonorM[1]))
        return 1;
    return 0;

}

bool Loop_For_Junctions(map <unsigned,unsigned char>::iterator &Start, map<unsigned,unsigned char>::iterator & End, unsigned char Sig, int Pairing, unsigned x,unsigned y, char* Read, Hash & Junctions, const int Direction)
{
    OPX* tempJuncList = new OPX[HIGHREPTHRES];
    unsigned lastY;
    unsigned tempJuncListPtr = 0;
    bool highRep = false;
    unsigned tempCounter = 0;
    bool Good_Junc = false;
    OP JPair;

    while(Start != End)
    {
        if (Start->second == Sig)
        {
            int L=0,Mis_Count=0;
            if(Direction == -4)
            {
                for (int i=y-4,j=Start->first-4;i && Mis_Count <= MIS_IN_UNMAPPED;i--,j--)
                {
                    if (Read[i] != Bit2Char(j)) Mis_Count++;
                    L++;
                }
            }
            else
            {
                for (int i=y+4,j=Start->first+4;i < STRINGLENGTH && Mis_Count <= MIS_IN_UNMAPPED;i++,j++)
                {
                    if (Read[i] != Bit2Char(j)) Mis_Count++;
                    L++;
                }
            }
            bool Junc_OK=false;
            bool addToList = true;
            if (L<= RESIDUE)
            {
                if (Mis_Count <=MIS_IN_RES) Junc_OK=true;
            }
            else if (Mis_Count<=MIS_IN_UNMAPPED) Junc_OK=true;
            if (Junc_OK)
                //if (Mis_Count<=MIS_IN_UNMAPPED)
            {
                //JPair.x=Start->first;
                //JPair.y=y;
                //Junctions.Insert(JPair,Pairing);
                //if (Pairing == S_ct_ac) {Good_Junc=true;break;}
                if(tempCounter == 0)
                {
					if(Direction == -4)
					{
						tempJuncList[0].x = Start->first;
						tempJuncList[0].y = x;
					}
					else
					{
                    	tempJuncList[0].x = x;
                    	tempJuncList[0].y = Start->first;
					}
					tempJuncList[0].Motif = Pairing;
                    tempJuncListPtr++;
                    lastY = Start->first;
                    addToList = false;
                }
                else if(tempCounter >= HIGHREPTHRES)
                {
                    printf("tempJuncList size too small\n");
                    highRep = true;
                    break;
                }
                else if(tempCounter > 0)
                {
                    if( dist(Start->first,lastY) < MAXMULTIGAP)
                    {
                        addToList = false;
                    }
                    if( tempJuncListPtr > 0 && dist(tempJuncList[tempJuncListPtr-1].y, Start->first) < MAXMULTIGAP)
                    {
                        tempJuncListPtr-=1;
                        addToList = false;
                    }
                }
                if(addToList)
                {
					if(Direction == -4)
					{
						tempJuncList[tempJuncListPtr].x = Start->first;
						tempJuncList[tempJuncListPtr].y = x;
					}
					else
					{
                    	tempJuncList[tempJuncListPtr].x = x;
                    	tempJuncList[tempJuncListPtr].y = Start->first;
					}
                    tempJuncList[tempJuncListPtr].Motif = Pairing;
                    tempJuncListPtr++;
                }
                lastY = Start->first;
                tempCounter++;
            }
        }
        Start++;
    }
    if(!highRep)
    {
        for(int i = 0; i<tempJuncListPtr;i++)
        {
            JPair.x = tempJuncList[i].x;
            JPair.y = tempJuncList[i].y;
            Junctions.Insert(JPair,tempJuncList[i].Motif);
        }
        if(tempJuncListPtr == 1)
            Good_Junc= true;
    }
    delete [] tempJuncList;
    return Good_Junc;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Match_Unmapped...
 *  Description:  Take unmapped read and try to map it.. Long -Short (+)[exon.bigright]
 * =====================================================================================
 */

void Match_Unmapped(FILE* Unmapped,Hash & Junctions,COVTYPE* Coverage,map <unsigned,unsigned char> & Acc_AG,map <unsigned,unsigned char> & Acc_AC)
{
	rewind(Unmapped);
	Unmapped_Record URecord;
	char Read[MAXTAG];
	int Mis_Count=0;
	unsigned char A,B,C,D,Sig;
	map<unsigned,unsigned char>::iterator Start,End,End_AG,End_AC; 
	End_AG=Acc_AG.end();
	End_AC=Acc_AC.end();

	while (1)
	{
		fread(&URecord,sizeof(URecord),1,Unmapped);
		STRINGLENGTH=URecord.Strlen;
		if(!fread(Read,STRINGLENGTH,1,Unmapped)) break;
		//for (unsigned x=URecord.Location+URecord.Length,y=URecord.Length+1;(x>URecord.Location+URecord.Length-EXTEND) && (y+RESIDUE<STRINGLENGTH);x--,y--)
		for (unsigned x=URecord.Location+URecord.Length,y=URecord.Length+1;(x>URecord.Location+URecord.Length-EXTEND) ;x--,y--)
		{
			if(y+RESIDUE>STRINGLENGTH) continue; //The residue is too small
			//if (Coverage[x] < COVERAGE_THRESH) continue;
			unsigned char DonorA= Bit2Char(x+1);
			unsigned char DonorB= Bit2Char(x+2);

			bool Good_Junc=false;
			int Pairing =0;
			if (DonorA == N_g)
			{
				if(DonorB == N_t)//GT [GT-AG]
				{
					A=Read[y];B=Read[y+1];C=Read[y+2];D=Read[y+3];
					Sig= (D | (C<<2) | (B<<4) |(A<<6));
					Start=Acc_AG.lower_bound(x);End=Acc_AG.upper_bound(x+EXONGAP);
                    Pairing = S_gt_ag;
                }
                else if(DonorB == N_c)//GC [GC-AG]
                {
					A=Read[y];B=Read[y+1];C=Read[y+2];D=Read[y+3];
					Sig= (D | (C<<2) | (B<<4) |(A<<6));
					Start=Acc_AG.lower_bound(x);End=Acc_AG.upper_bound(x+EXONGAP);
                    Pairing = S_gc_ag;
                }
			    if (End_AG==Start) Start=End;//End is the position after the last element;
            }
            else if(DonorA == N_a && DonorB == N_t)//AT [AT-AC]
            {
				A=Read[y];B=Read[y+1];C=Read[y+2];D=Read[y+3];
				Sig= (D | (C<<2) | (B<<4) |(A<<6));
				Start=Acc_AC.lower_bound(x);End=Acc_AC.upper_bound(x+EXONGAP);
                Pairing = S_at_ac;
                if (End_AC == Start) Start = End;
            }
            if(Pairing)
                Good_Junc = Loop_For_Junctions(Start, End, Sig,Pairing,x,y,Read,Junctions,4);
			if (Good_Junc) break;
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Match_Unmapped...
 *  Description:  Take unmapped read and try to map it.. Short Left-Long Right (-) [exonm.bigright]
 * =====================================================================================
 */
void Match_UnmappedMX(FILE* Unmapped,Hash & Junctions,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Donor_CT)
{
    rewind(Unmapped);
    Unmapped_Record URecord;
    char Read[MAXTAG];
    int Mis_Count=0;
    unsigned char A,B,C,D,Sig;
    map<unsigned,unsigned char>::iterator Start,End,End_CT,End_GT; 
    End_CT=Donor_CT.end(),End_GT=Donor_GT.end();

    while (1)
    {
        fread(&URecord,sizeof(URecord),1,Unmapped);
        STRINGLENGTH=URecord.Strlen;
        if(!fread(Read,STRINGLENGTH,1,Unmapped)) break;
        //for (unsigned y=URecord.Location,x=STRINGLENGTH-URecord.Length-1;(y<URecord.Location+EXTEND) && x>RESIDUE-1;x++,y++)
        for (unsigned y=URecord.Location,x=STRINGLENGTH-URecord.Length-1;(y<URecord.Location+EXTEND) ;x++,y++)
        {
            if(x<RESIDUE-1) continue;
            unsigned char AccA= Bit2Char(y-2);
            unsigned char AccB= Bit2Char(y-1);

            bool Good_Junc=false;
            int Pairing =0;
            if (AccB == N_c)
            {
                if(AccA == N_a||AccA == N_g)// AC/GC [GC-AG]-->[CT-GC],[GT-AG]-->[CT - AC]
                {
                    if (AccA == N_a) Pairing=S_ct_ac; else Pairing=S_ct_gc;
                    A=Read[x-3];B=Read[x-2];C=Read[x-1];D=Read[x];
                    Sig= (D | (C<<2) | (B<<4) |(A<<6));
                    Start=Donor_CT.lower_bound(y-EXONGAP);End=Donor_CT.upper_bound(y);
                    if (End_CT==Start) Start=End;
                }
            }
            else if(AccA==N_a && AccB == N_t)//GT [AT-AC]-->[GT - AT]
            {
                A=Read[x-3];B=Read[x-2];C=Read[x-1];D=Read[x];
                Sig= (D | (C<<2) | (B<<4) |(A<<6));
                Start=Donor_GT.lower_bound(y-EXONGAP);End=Donor_GT.upper_bound(y);
                Pairing = S_gt_at;
                if (End_GT==Start) Start=End;
            }
            if(Pairing)
                Good_Junc = Loop_For_Junctions(Start,End,Sig,Pairing,y,x,Read,Junctions,-4);
            if(Good_Junc) break;
        }
    }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Match_Unmapped...
 *  Description:  Take unmapped read and try to map it..  Short Left - Long Right (+)[exon.bigright]
 * =====================================================================================
 */
void Match_UnmappedX(FILE* Unmapped,Hash & Junctions,map <unsigned,unsigned char> & Donor_AT,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Donor_GC)
{
	rewind(Unmapped);
	Unmapped_Record URecord;
	char Read[MAXTAG];
	int Mis_Count=0;
	unsigned char A,B,C,D,Sig;
	map<unsigned,unsigned char>::iterator Start,End,End_GT,End_GC,End_AT; 
	End_GT=Donor_GT.end();
	End_GC=Donor_GC.end();
	End_AT=Donor_AT.end();

	while (1)
	{
		fread(&URecord,sizeof(URecord),1,Unmapped);
		STRINGLENGTH=URecord.Strlen;
		if(!fread(Read,STRINGLENGTH,1,Unmapped)) break;
		//for (unsigned y=URecord.Location,x=STRINGLENGTH-URecord.Length-1;(y<URecord.Location+EXTEND) && x>RESIDUE -1;x++,y++)
		for (unsigned y=URecord.Location,x=STRINGLENGTH-URecord.Length-1;(y<URecord.Location+EXTEND);x++,y++)
		{
			if (x<RESIDUE-1) continue;
			unsigned char AccA= Bit2Char(y-2);
			unsigned char AccB= Bit2Char(y-1);

			int Pairing =0;
			bool Good_Junc=false;
			if (AccA == N_a)
			{
				if(AccB == N_g)//AG [GT/GC-AG]
				{
					A=Read[x-3];B=Read[x-2];C=Read[x-1];D=Read[x];
                    Sig= (D | (C<<2) | (B<<4) |(A<<6));
                    Start=Donor_GT.lower_bound(y-EXONGAP);End=Donor_GT.upper_bound(y);
                    if (End_GT==Start) Start=End;
                    Pairing = S_gt_ag;
                    Good_Junc = Loop_For_Junctions(Start,End,Sig,Pairing,y,x,Read,Junctions,-4);
                    if(Good_Junc) {break;}

                    Start=Donor_GC.lower_bound(y-EXONGAP);End=Donor_GC.upper_bound(y);
                    Pairing = S_gc_ag;
                    if (End_GC==Start) Start=End;
                    Good_Junc = Loop_For_Junctions(Start,End,Sig,Pairing,y,x,Read,Junctions,-4);
                    if(Good_Junc)
                    {
                        break;
                    }
                }
                
				else if(AccB == N_c)//AT [AC-AT]
				{
					A=Read[x-3];B=Read[x-2];C=Read[x-1];D=Read[x];
					Sig= (D | (C<<2) | (B<<4) |(A<<6));
					Start=Donor_AT.lower_bound(y-EXONGAP);End=Donor_AT.upper_bound(y);
                    Pairing = S_at_ac;
                    if (End_AT==Start) Start=End;
                    Good_Junc = Loop_For_Junctions(Start,End,Sig,Pairing,y,x,Read,Junctions,-4);
                    if(Good_Junc)
                    {
                        break;
                    }
                }
			} 
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Match_Unmapped...
 *  Description:  Take unmapped read and try to map it.. Long -Short (-)
 * =====================================================================================
 */
void Match_UnmappedM(FILE* Unmapped,Hash & Junctions,map <unsigned,unsigned char> & Acc_GC,map <unsigned,unsigned char> & Acc_AC,map <unsigned,unsigned char> & Acc_AT)
{
	rewind(Unmapped);
	Unmapped_Record URecord;
	char Read[MAXTAG];
	int Mis_Count=0;
	static unsigned Debug_Count=0;
	unsigned char A,B,C,D,Sig;
	map<unsigned,unsigned char>::iterator Start,End,End_AAC,End_AGC,End_AAT; 
	End_AAC=Acc_AC.end();
	End_AGC=Acc_GC.end();
	End_AAT=Acc_AT.end();

	while (1)
	{
		fread(&URecord,sizeof(URecord),1,Unmapped);
		STRINGLENGTH=URecord.Strlen;
		if(!fread(Read,STRINGLENGTH,1,Unmapped)) break;
//Donor in complement..
		//for (unsigned x=URecord.Location+URecord.Length,y=URecord.Length+1;(x>URecord.Location+URecord.Length-EXTEND) && (y+RESIDUE<STRINGLENGTH);x--,y--)
		for (unsigned x=URecord.Location+URecord.Length,y=URecord.Length+1;(x>URecord.Location+URecord.Length-EXTEND);x--,y--)
		{
			if(y+RESIDUE>STRINGLENGTH) continue;
			unsigned char DonorA= Bit2Char(x+1);//(unsigned char)(*(Original_Text+(Orig_Text+x+1)/4) << (((Orig_Text+x+1) % 4) * 2))>>6;
			unsigned char DonorB= Bit2Char(x+2);//(unsigned char)(*(Original_Text+(Orig_Text+x+2)/4) << (((Orig_Text+x+2) % 4) * 2))>>6;

			int Pairing =0;
			bool Good_Junc=false;
			if (DonorB == N_t)
			{
				if(DonorA == N_c)//CT [GT-AG]-->CT - AC
				{
					A=Read[y];B=Read[y+1];C=Read[y+2];D=Read[y+3];
					Sig= (D | (C<<2) | (B<<4) |(A<<6));
					Start=Acc_AC.lower_bound(x);End=Acc_AC.upper_bound(x+EXONGAP);
                    Pairing = S_ct_ac;
					if (End_AAC==Start) Start=End;
                    Good_Junc = Loop_For_Junctions(Start,End,Sig,Pairing,x,y,Read,Junctions,4);
                    if (Good_Junc){ break;}

					Start=Acc_GC.lower_bound(x);End=Acc_GC.upper_bound(x+EXONGAP);
                    Pairing = S_ct_gc;
					if (End_AGC==Start) Start=End;
                    Good_Junc = Loop_For_Junctions(Start,End,Sig,Pairing,x,y,Read,Junctions,4);
                    if(Good_Junc) {break;}
				}
				else if(DonorA == N_g)//gt [at-ac]--> gt-at
				{
					A=Read[y];B=Read[y+1];C=Read[y+2];D=Read[y+3];
					Sig= (D | (C<<2) | (B<<4) |(A<<6));
					Start=Acc_AT.lower_bound(x);End=Acc_AT.upper_bound(x+EXONGAP);
                    Pairing = S_gt_at;
					if (End_AAT==Start) Start=End;
					Good_Junc = Loop_For_Junctions(Start,End,Sig,Pairing,x,y,Read,Junctions,4);
                    if(Good_Junc) {break; }
				}
			} 
		}
	}

}/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Seek_Island_Junctions...
 *  Description:  Extend islands and seek splice signals..
 * =====================================================================================
 */
void Seek_Island_Junctions(Island* Island_List,map <unsigned,unsigned char> & Donor_AT,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Donor_GC,map <unsigned,unsigned char> & Acc_AG,map <unsigned,unsigned char> & Acc_AC)
{
//[GT-AG],[GC - AG],[AC-AT]

	Island* Current_Island=Island_List;
	Island* Next_Island;
	while (Current_Island)
	{
		if (Current_Island->Status != DELETED)
		{
			//------------------------------ Seek canonical Acceptor sites ------------------------------------------------------
			for (unsigned y=Current_Island->Start+EXTENDISLAND;y>Current_Island->Start-EXTENDISLAND;y--)
			//for (unsigned y=Current_Island->End+EXTENDISLAND;y>Current_Island->Start-EXTENDISLAND;y--)
			{

				unsigned char AccA= Bit2Char(y-2);//(unsigned char)(*(Original_Text+(Orig_Text+y-2)/4) << (((Orig_Text+y-2) % 4) * 2))>>6;
				unsigned char AccB= Bit2Char(y-1);//(unsigned char)(*(Original_Text+(Orig_Text+y-1)/4) << (((Orig_Text+y-1) % 4) * 2))>>6;


				int Pairing =0;
				if (AccA == N_a)
				{
					if(AccB == N_g)//AG -high chance case [GT-AG] or [GC - AG]
					{ 
						unsigned char A,B,C,D;
						A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
						Acc_AG[y]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing = S_ag;
					}
					else if (AccB == N_c)//AC - [AT-AC]
					{
						unsigned char A,B,C,D;
						A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
						Acc_AC[y]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing = S_ac;
					}
				}
				//if (Pairing) printf("[%d\n",y);
			}
			//------------------------------ Seek canonical Donor sites ------------------------------------------------------
			for (unsigned x=Current_Island->End-EXTENDISLAND;x<Current_Island->End+EXTENDISLAND;x++)
			//for (unsigned x=Current_Island->Start-EXTENDISLAND;x<Current_Island->End+EXTENDISLAND;x++)
			{
				unsigned char DonorA= Bit2Char(x+1);//(unsigned char)(*(Original_Text+(Orig_Text+x+1)/4) << (((Orig_Text+x+1) % 4) * 2))>>6;
				unsigned char DonorB= Bit2Char(x+2);//(unsigned char)(*(Original_Text+(Orig_Text+x+2)/4) << (((Orig_Text+x+2) % 4) * 2))>>6;

				int Pairing =0;
				if (DonorA == N_g)
				{
					if(DonorB == N_t)//GT [GT-AG]
					{
						unsigned char A,B,C,D;
						A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
						Donor_GT[x]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing= S_gt;
					}
					else if(DonorB == N_c)//GC [GC-AG]
					{
						unsigned char A,B,C,D;
						A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
						Donor_GC[x]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing= S_gc;
					}
				} 
				if (DonorA == N_a && DonorB == N_t)//AT [AT-AC]
				{
					unsigned char A,B,C,D;
					A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
					Donor_AT[x]= (D | (C<<2) | (B<<4) |(A<<6));
					Pairing = S_at; 
				}
				//if (Pairing) printf("%d]\n",x);

			}
		}
		Current_Island=Current_Island->Next;
	}
}

void Seek_Island_JunctionsM(Island* Island_List,map <unsigned,unsigned char> & Donor_CT,map <unsigned,unsigned char> & Donor_GT,map <unsigned,unsigned char> & Acc_GC,map <unsigned,unsigned char> & Acc_AC,map <unsigned,unsigned char> & Acc_AT)
{
//[CT-AC],[CT - GC],[AT-GT]

	Island* Current_Island=Island_List;
	Island* Next_Island;
	while (Current_Island)
	{
		if (Current_Island->Status != DELETED)
		{
			//------------------------------ Seek canonical Acceptor sites ------------------------------------------------------
			for (unsigned y=Current_Island->Start+EXTENDISLAND;y>Current_Island->Start-EXTENDISLAND;y--)
			//for (unsigned y=Current_Island->End+EXTENDISLAND;y>Current_Island->Start-EXTENDISLAND;y--)
			{

				unsigned char AccA= Bit2Char(y-2);//(unsigned char)(*(Original_Text+(Orig_Text+y-2)/4) << (((Orig_Text+y-2) % 4) * 2))>>6;
				unsigned char AccB= Bit2Char(y-1);//(unsigned char)(*(Original_Text+(Orig_Text+y-1)/4) << (((Orig_Text+y-1) % 4) * 2))>>6;


				int Pairing =0;
				if (AccA == N_a)
				{
					if(AccB == N_c)//GC [CT-GC]
					{ 
						unsigned char A,B,C,D;
						A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
						Acc_AC[y]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing = S_ac;
					}
					else if (AccB == N_t)//GT - [AT-GT]
					{
						unsigned char A,B,C,D;
						A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
						Acc_AT[y]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing = S_at;
					}
				}
				if (AccA == N_g && AccB == N_c)//AC [CT-AC]
				{
					unsigned char A,B,C,D;
					A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
					Acc_GC[y]= (D | (C<<2) | (B<<4) |(A<<6));
					Pairing = S_gc; 
				}
				//if (Pairing) printf("[%d\n",y);
			}
			//------------------------------ Seek canonical Donor sites ------------------------------------------------------
			for (unsigned x=Current_Island->End-EXTEND;x<Current_Island->End+EXTEND;x++)
			//for (unsigned x=Current_Island->Start-EXTENDISLAND;x<Current_Island->End+EXTENDISLAND;x++)
			{
				unsigned char DonorA= Bit2Char(x+1);//(unsigned char)(*(Original_Text+(Orig_Text+x+1)/4) << (((Orig_Text+x+1) % 4) * 2))>>6;
				unsigned char DonorB= Bit2Char(x+2);//(unsigned char)(*(Original_Text+(Orig_Text+x+2)/4) << (((Orig_Text+x+2) % 4) * 2))>>6;

				int Pairing =0;
				if (DonorB == N_t)
				{
					if(DonorA == N_c)//CT [CT-AC][CT -GC]
					{
						unsigned char A,B,C,D;
						A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
						Donor_CT[x]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing= S_ct;
					}
					else if(DonorA == N_g)//AT [AT-GT]
					{
						unsigned char A,B,C,D;
						A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
						Donor_GT[x]= (D | (C<<2) | (B<<4) |(A<<6));
						Pairing= S_gt;
					}
				} 
				//if (Pairing) printf("%d]\n",x);

			}
		}
		Current_Island=Current_Island->Next;
	}
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Print_BED
 *  Description:  Print the coverage wiggle...
 * =====================================================================================
 */
void Print_BED(Hash & Junctions,Island* Island_List, COVTYPE* Coverage,char* chromosome,char Strand)
{
	static char NOHEADER=TRUE;
	printf("\t+Writing Files..\n");
//-------------------------------------- Wigglegram ------------------------------------------------------------------------------------------


	Island* Current_Island=Island_List;
	Island* Next_Island;
	while (Current_Island->Status== DELETED){if(!(Current_Island=Current_Island->Next)) break;}//skip deleted islands...
	while (Current_Island)
	{
		unsigned i=Current_Island->Start;
		for (;i<Current_Island->End ;)
		{
			unsigned First_Base=i;
			while(i<Current_Island->End && Coverage[First_Base]==Coverage[++i]);
			fprintf (BigBed,"%s\t%u\t%u\t%u\n",chromosome,First_Base,i,Coverage[First_Base]);
		}

		if(Current_Island=Current_Island->Next)
		{
			while (Current_Island->Status== DELETED){if(!(Current_Island=Current_Island->Next)) break;}//skip deleted islands...
			if(Current_Island) fprintf (BigBed,"%s\t%u\t%u\t%u\n",chromosome,i,Current_Island->Start,0);// else fprintf (BigBed,"%s\t%u\t%u\t%u\n",chromosome,i,Current_Island->Start,0);
		}
	}

//-------------------------------------- Wigglegram ------------------------------------------------------------------------------------------
//-------------------------------------- Junctions ------------------------------------------------------------------------------------------

	static int j=0;
	JStat JStat;
	OP JPair;
	map <unsigned,unsigned>::iterator ExonH_I;
	map <unsigned,unsigned>::iterator ExonT_I;

	if (NOHEADER) {fprintf(Junctions_File,"track name=\"Rawbin Exons\"\n");NOHEADER=FALSE;}
	char Junc_Not_Empty = Junctions.Init_Iterate(JPair,JStat);

	while (Junc_Not_Empty)
	{
		if(JStat.Junc_Type)
		{

			OP H,T;

			char Strand_Sign = (JStat.Junc_Type >= MINUS_JUNC) ? '-' : '+';
			int Left_Cover=Coverage[JPair.x];int Right_Cover=Coverage[JPair.y];
			int Left_Drop=Coverage[JPair.x]-Coverage[JPair.x+1];int Right_Drop=Coverage[JPair.y]-Coverage[JPair.y-1];
			H.x=JPair.x-2;H.y=JPair.x+1;
			T.x=JPair.y;T.y=JPair.y+2;
			if(H.x>=T.x || T.x >= T.y || H.x >= H.y)
			{
				printf ("Print_BED():Junction Parse error..\nH.x:H.y=[%u,%u]       T.x:T.y=[%u,%u]\n",H.x,H.y,T.x,T.y);exit(0);
			}
			bool Pass_Junc=true;
			if(H.y+EX_MIN >= T.x) Pass_Junc=false;//Junctions are not too close.. 
			if (JStat.Junc_Type != S_gt_ag && JStat.Junc_Type != S_ct_ac)//require strong evidence for rare junctions..
			{
				if (Left_Cover < COVERAGE_THRESH || Right_Cover < COVERAGE_THRESH) Pass_Junc=false;//Junction is loosely covered..
				if (Left_Drop <=0 || Right_Drop<=0) Pass_Junc=false;//junction does not have peak coverage..
			}
			/*   
			     1.  chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
			     2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
			     3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
			     4. name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
			     5. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
			     shade 	  	  	  	  	  	  	  	  	 
			     6. strand - Defines the strand - either '+' or '-'.
			     7. thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
			     8. thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
			     9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
			     10. blockCount - The number of blocks (exons) in the BED line.
			     11. blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
			     12. blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 
			*/
			//if ((JStat.Junc_Type-255 !=1) && (JStat.Junc_Type-255 !=4)) Pass_Junc=false;
			if(Pass_Junc)
			{
				char const *Junction_Type=Get_Junc_Type(JStat.Junc_Type);
				fprintf (Junctions_File,\
						"%s\t%u\t%u\tJUNC%d_%s_%d_%d_%d_%d\t%u\t%c\t%u\t%u\t255,0,0\t2\t%u,%u\t%u,%u\n",
						chromosome,
						H.x,T.y+1, //Chrom start and end...
						j++,//Junction name..
						Junction_Type,Left_Cover,Right_Cover,Left_Drop,Right_Drop,
						JStat.Count,//score to color...
						Strand_Sign,
						H.x,T.y,//thickStart and end
						//H.y-H.x+1,T.y-T.x+1,//block sizes..
						H.y-H.x,T.y-T.x+1,//block sizes..
						0,T.x-H.x//block starts...
					);
			}

		}
		Junc_Not_Empty=Junctions.Iterate(JPair,JStat);
	}
}

inline char* Get_Junc_Type(int Junc_Type)
{
	return (char*) Int_To_Junc[Junc_Type-255];
}
//-------------------------------------- Exons ------------------------------------------------------------------------------------------
//}--------------------------------------  Coverage -------------------------------------------------------------------

//{--------------------------------------  Two Mismatch -------------------------------------------------------------------

void Search_Backwards_X10_OneSA(struct SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int & Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}


	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fwfmi->inverseSa0) Index--;//adjust for missing $
		Now=fwfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fwfmi->cumulativeFreq[Now] + BWTOccValue(fwfmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)//only one mismatch allowed here...
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					if(!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=1;
					Search_Backwards(Tag,2,LH,LH,Two_Mismatches_At_End_Pointer,Two_Mismatches_At_End,Mismatches_Backward_Pointer,Mismatches_Backward);
					return;
				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else return;
	}
}

void Search_Backwards_X10(struct SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int & Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward)
{
	if (!Tag.Start) return;
	int BMHStack_Top=0;
	BMHStack[0]=Tag;
	struct SARANGE Range,Temp_Range;
	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_X10_OneSA(Range,Count,Start,StringLength,Two_Mismatches_At_End_Pointer,Two_Mismatches_At_End,Mismatches_Backward_Pointer,Mismatches_Backward);
			if(MAXHITS==HITS) return;
		}
		else
		{
			Branch_Detect_Backwards(Range,fwfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								Temp_Range.Level=1;
								Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
								memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARANGE));
								Search_Backwards(Temp_Range,2,LH,LH,Two_Mismatches_At_End_Pointer,Two_Mismatches_At_End,Mismatches_Backward_Pointer,Mismatches_Backward);
								if(MAXHITS==HITS) return;
								Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
								memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARANGE));
							}
							else continue;
						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
				} 
			}
		}
	}
	return;
}

void Reverse(struct SARANGE & Tag,int Start,int StringLength)
{	
	unsigned Temp=0;
	char New_Char;
	char Mismatch_Count=Tag.Mismatches;
	unsigned pos;
	for( int i=0;i<Mismatch_Count;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Temp=Temp | (Current_Tag[pos]<<i*2);
		Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
	}
	Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
	{
		Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;
		Search_Backwards_Exact(Tag,Current_Tag,STRINGLENGTH,RH);
	}
	Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
	for( int i=0;i<Tag.Mismatches;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Current_Tag[pos]=(Temp>>(2*i)) & 3;
	}
	return;
}

void Search_X01_OneSA(SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int & Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= revfmi->inverseSa0) Index--;//adjust for missing $
		Now=revfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					Reverse(Tag,STRINGLENGTH,RH);
					if(Tag.Start)
					{
						Tag.Level=1;
						Search_Backwards(Tag,2,LH,LH,Two_Mismatches_At_End_Pointer,Two_Mismatches_At_End,Mismatches_Backward_Pointer,Mismatches_Backward);
					}
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else return;
	}
}

void Search_X01(SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int & Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	BMHStack[0]=Tag;
	struct SARANGE Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_X01_OneSA(Range,Count,Start,StringLength,Two_Mismatches_At_End_Pointer,Two_Mismatches_At_End,Mismatches_Backward_Pointer,Mismatches_Backward);
			if(MAXHITS==HITS) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								Reverse(Temp_Range,STRINGLENGTH,RH);
								if(Temp_Range.Start)
								{
									Temp_Range.Level=1;
									Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
									Search_Backwards(Temp_Range,2,LH,LH,Two_Mismatches_At_End_Pointer,Two_Mismatches_At_End,Mismatches_Backward_Pointer,Mismatches_Backward);
									if(MAXHITS==HITS) return;
									Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
								}
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							BMHStack[FSHStack_Top]=Temp_Range;
						}
					}
				} 
			}

		}
	}
	return;
}


char Two_Mismatch_Search()
{
	if(Two_Mismatches_At_End_Forward_PointerP)//give priority to forward direction as most erros occur in the end..
	{
		Current_Tag=Head.Tag;
		Two_Mismatches_At_End_ForwardP[0].Mismatch_Pos[1]= (STRINGLENGTH-1);//mismatches of the form 0|2, with last mismatch at the end...
		Print_LocationX(Two_Mismatches_At_End_ForwardP[0],FW);
		HITS++;
		//if(MAXHITS==Hits) return TRUE;
		return TRUE;
	}
	if(Two_Mismatches_At_End_Forward_PointerC)//give priority to forward direction as most erros occur in the end..
	{
		Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH);
		Two_Mismatches_At_End_ForwardC[0].Mismatch_Pos[1]= (STRINGLENGTH-1);//mismatches of the form 0|2, with last mismatch at the end...
		Print_LocationX(Two_Mismatches_At_End_ForwardC[0],FW);
		HITS++;
		//if(MAXHITS==Hits) return TRUE;
		return TRUE;
	}

	if(Two_Mismatches_At_End_PointerP)
	{
		Current_Tag=Head.Tag;
		Print_LocationX(Two_Mismatches_At_EndP[0],BW);//Mismatches of the form 2|0, with one mismatch at the first position
		HITS++;
		//if(MAXHITS==Hits) return TRUE;
		return TRUE;
	}
	if(Two_Mismatches_At_End_PointerC)
	{
		Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH);
		Print_LocationX(Two_Mismatches_At_EndC[0],BW);//Mismatches of the form 2|0, with one mismatch at the first position
		HITS++;
		//if(MAXHITS==Hits) return TRUE;
		return TRUE;
	}

	if(Mismatches_Forward_PointerP)
	{
		Current_Tag=Head.Tag;
		for(int i=Mismatches_Forward_PointerP-1;i>=0;i--)
		{
			if(Search_Forwards(Mismatches_ForwardP[i],2,LH+1,RH,Two_Mismatches_At_End_Forward_PointerP,Two_Mismatches_At_End_ForwardP,Mismatches_Forward_PointerP,Mismatches_ForwardP)) return TRUE;
		}
	}

	if(Mismatches_Forward_PointerC)
	{
		Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH);
		for(int i=Mismatches_Forward_PointerC-1;i>=0;i--)
		{
			if(Search_Forwards(Mismatches_ForwardC[i],2,LH+1,RH,Two_Mismatches_At_End_Forward_PointerC,Two_Mismatches_At_End_ForwardC,Mismatches_Forward_PointerC,Mismatches_ForwardC)) return TRUE;
		}
	}

	if(Mismatches_Backward_PointerP)
	{
		Current_Tag=Head.Tag;
		for(int i=Mismatches_Backward_PointerP-1;i>=0;i--)
		{
			if (Search_Backwards(Mismatches_BackwardP[i],2,LH,LH,Two_Mismatches_At_End_PointerP,Two_Mismatches_At_EndP,Mismatches_Backward_PointerP,Mismatches_BackwardP)) return TRUE;
		}
	}

	if(Mismatches_Backward_PointerC)
	{
		Current_Tag=Head.Complement+(STRINGLENGTHO-STRINGLENGTH);
		for(int i=Mismatches_Backward_PointerC-1;i>=0;i--)
		{
			if (Search_Backwards(Mismatches_BackwardC[i],2,LH,LH,Two_Mismatches_At_End_PointerC,Two_Mismatches_At_EndC,Mismatches_Backward_PointerC,Mismatches_BackwardC)) return TRUE;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------------------
	for (int i=0;i<2;i++)
	{
		SARANGE Tag;
		Tag.Start=1;Tag.Level=1;Tag.End=SOURCELENGTH;Tag.Skip=0;Tag.Mismatches=0;Tag.Mismatch_Char=0;Current_Tag=Orientation[i].String;
		Search_Backwards_Exact(Tag,Current_Tag,STRINGLENGTH,RHQR);
		Tag.Level=1;

		Search_Backwards_X10(Tag,1,LH + RHQL, RHQL,Two_Mismatches_At_End_PointerP,Two_Mismatches_At_EndP,Mismatches_Backward_PointerP,Mismatches_BackwardP);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==HITS) return TRUE;

		Tag.Start=1;Tag.Level=1;Tag.End=SOURCELENGTH;Tag.Skip=0;Tag.Mismatches=0;Tag.Mismatch_Char=0;Current_Tag=Orientation[i].String;
		Search_Forwards_Exact(Tag,Current_Tag,LH+1,RHQL);
		Tag.Level=1;
		Search_X01(Tag,1,LH + RHQL +1,RHQR,Two_Mismatches_At_End_PointerP,Two_Mismatches_At_EndP,Mismatches_Backward_PointerP,Mismatches_BackwardP);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==HITS) return TRUE;
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
	return FALSE;
}

//}--------------------------------------  Two Mismatch -------------------------------------------------------------------

//{-------------------------------------- Mismatch Search Forward -------------------------------------------------------------------
void Branch_Detect (const struct SARANGE Tag,BWT *fmi,int Start)
{

	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;
	if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
	{
		unsigned Last, First;
		char Now;

		if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 

		for (unsigned long Pos=First;Pos<=Last;Pos++)
		{
			Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
			Branch_Characters[Now]++;	
		}

		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Tag.Level+Start] && Branch != Current_Tag[Start+Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else*/ 
			if (Branch_Characters[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Branch_Characters[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Tag.Level+Start] && Branch != Current_Tag[Start+Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else*/
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.End+1, Branch);
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
			}
		}

	}
}

char Search_Forwards(SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Forward_Pointer,SARANGE* Two_Mismatches_At_End_Forward,int & Mismatches_Forward_Pointer,SARANGE* Mismatches_Forward)
//char Search_Backwards(SARANGE & Tag,int Count,int Start,int StringLength,int Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward);
{
	if (!Tag.Start) return FALSE;
	Start=Start-2;//Adjust for offset difference
	int FSStack_Top=0;
	FSSStack[0]=Tag;
	struct SARANGE Range,Temp_Range;
	while(FSStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSSStack[FSStack_Top];
		FSStack_Top--;		//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Forwards_OneSA(Range,Count,Start,StringLength,Two_Mismatches_At_End_Forward_Pointer,Two_Mismatches_At_End_Forward,Mismatches_Forward_Pointer,Mismatches_Forward);
			if(MAXHITS==HITS) return TRUE;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{

						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;

					}


					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (Temp_Range.Mismatches == Count) //dont print exact matches
							{
								Print_LocationX(Temp_Range,FW);
								HITS++;
								//if (MAXHITS==HITS) return TRUE;
								return TRUE;
							}
							else continue;
						}
						else 
						{

							FSStack_Top++;//Push range
							Temp_Range.Level++;
							FSSStack[FSStack_Top]=Temp_Range;
						}
					}
					else
					{
						if(2 >Count)//store only for one mismatch... and last node will not branch
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							else //2 mismatches with the last at the end...
							{
								if(Two_Mismatches_At_End_Forward_Pointer < END_BOUND)
								{
									Two_Mismatches_At_End_Forward[Two_Mismatches_At_End_Forward_Pointer]=Temp_Range;
									Two_Mismatches_At_End_Forward_Pointer++;
								}
								continue;
							}
							if(Mismatches_Forward_Pointer < ARRAY_BOUND)
							{
								Mismatches_Forward[Mismatches_Forward_Pointer]=Temp_Range;
								Mismatches_Forward_Pointer++;
							}
						}
						continue;
					}
				} 
			}
		}
	}
	return FALSE;
}

//void Search_Forwards_OneSA(struct SARANGE & Tag,int Count,int Start,int StringLength)
void Search_Forwards_OneSA(SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Forward_Pointer,SARANGE* Two_Mismatches_At_End_Forward,int & Mismatches_Forward_Pointer,SARANGE* Mismatches_Forward)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		
		Index=Tag.Start;
		if (Index >= revfmi->inverseSa0) Index--;//adjust for missing $
		Now=revfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		//if (!Do_Branch[Tag.Level+Start] && Current_Tag[Tag.Level+Start]!=Now) return;  
		Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Tag.Level+Start]!=Now)
		{

			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
			
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)//avoid printing exact matches...
				{
					//if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					if (!Tag.Skip) Tag.End=Tag.Start;
					Print_LocationX(Tag,FW);
					HITS++;
				}
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else//log 2 mismatches 
		{
			if(2 > Count)//store only for one mismatch... 
			{
				if (!Tag.Skip) Tag.End=Tag.Start;//possibly two mismatch exists..
				if (Tag.Level != StringLength) Tag.Level++; 
				else //2 mismatches occuring in last position...
				{
					if(Two_Mismatches_At_End_Forward_Pointer < END_BOUND)
					{
						//if(Tag.Skip) Tag.Start=Tag.End;
						Two_Mismatches_At_End_Forward[Two_Mismatches_At_End_Forward_Pointer]=Tag;
						Two_Mismatches_At_End_Forward_Pointer++;
					}
					return;
				}
				if(Mismatches_Forward_Pointer < ARRAY_BOUND)
				{
					Mismatches_Forward[Mismatches_Forward_Pointer]=Tag;
					Mismatches_Forward_Pointer++;
				}
			}
			return;
		}
	}
}
//}-------------------------------------- Mismatch Search -------------------------------------------------------------------

//{-------------------------------------- Mismatch Search backward -------------------------------------------------------------------

void Branch_Detect_Backwards (const struct SARANGE Tag,BWT *fmi,int Start)
{

	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;
	if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
	{
		unsigned Last, First;
		char Now;

		if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 

		for (unsigned long Pos=First;Pos<=Last;Pos++)
		{
			Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
			Branch_Characters[Now]++;	
		}

		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Start-Tag.Level] && Branch != Current_Tag[Start-Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else*/
			if (Branch_Characters[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Branch_Characters[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Start-Tag.Level] && Branch != Current_Tag[Start-Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else*/
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.End+1, Branch);
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
			}
		}

	}
}

//void Search_Backwards_OneSA(struct SARANGE & Tag,int Count,int Start,int StringLength)
void Search_Backwards_OneSA(SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int & Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fwfmi->inverseSa0) Index--;//adjust for missing $
		Now=fwfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if (!Do_Branch[Start-Tag.Level] && Current_Tag[Start-Tag.Level]!=Now) return;  
		Tag.Start = fwfmi->cumulativeFreq[Now] + BWTOccValue(fwfmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)
				{
					//if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					if (!Tag.Skip) Tag.End=Tag.Start;
					Print_LocationX(Tag,BW);
					HITS++;
				}
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else 
		{
			if(5 >= Tag.Mismatches && 5 > Count)// 2 mismatches
			{
				if(!Tag.Skip) Tag.End=Tag.Start;//possibly two mismatch exists..
				if (Tag.Level != StringLength) Tag.Level++; 
				else//two mismatches with the last at the end ... 
				{
					//if(Tag.Skip) Tag.Start=Tag.End;
					if(Two_Mismatches_At_End_Pointer < END_BOUND)
					{
						Two_Mismatches_At_End[Two_Mismatches_At_End_Pointer]=Tag;
						Two_Mismatches_At_End_Pointer++;
					}
					return;
				}
				if(Mismatches_Backward_Pointer < ARRAY_BOUND)
				{
					Mismatches_Backward[Mismatches_Backward_Pointer]=Tag;
					Mismatches_Backward_Pointer++;
				}
			}
			return;
		} 
	}
}

char Search_Backwards(struct SARANGE & Tag,int Count,int Start,int StringLength,int & Two_Mismatches_At_End_Pointer,SARANGE* Two_Mismatches_At_End,int & Mismatches_Backward_Pointer,SARANGE* Mismatches_Backward)
{
	if (!Tag.Start) return FALSE;
	int BMStack_Top=0;
	BMStack[0]=Tag;
	struct SARANGE Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMStack[BMStack_Top];
		BMStack_Top--;	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_OneSA(Range,Count,Start,StringLength,Two_Mismatches_At_End_Pointer,Two_Mismatches_At_End,Mismatches_Backward_Pointer,Mismatches_Backward);
			if(MAXHITS==HITS) return TRUE;
		}
		else
		{
			Branch_Detect_Backwards(Range,fwfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=Start-Temp_Range.Level;
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches==Count)
							{
								Print_LocationX(Temp_Range,BW);
								HITS++;
								//if(MAXHITS==HITS) return TRUE;
								return TRUE;
							}
							else continue;
						}
						else
						{
							BMStack_Top++;//Push range
							Temp_Range.Level++;
							BMStack[BMStack_Top]=Temp_Range;
						}
					}
					else 
					{
						if(2 > Count)// 2 mismatches...
						{
							if(Temp_Range.Level != StringLength) Temp_Range.Level++;
							//if((Start-Temp_Range.Level) != 0) Temp_Range.Level++;
							else // 2 mismatches with the last at the end?
							{
								if(Two_Mismatches_At_End_Pointer < END_BOUND)
								{
									Two_Mismatches_At_End[Two_Mismatches_At_End_Pointer]=Temp_Range;
									Two_Mismatches_At_End_Pointer++;
								}
								continue;
							}
							if (Mismatches_Backward_Pointer < ARRAY_BOUND)
							{
								Mismatches_Backward[Mismatches_Backward_Pointer]=Temp_Range;
								Mismatches_Backward_Pointer++;
							}
						}
						continue;
					}
				} 
			}
		}
	}
	return FALSE; 
}

//}-------------------------------------- Mismatch Search backward -------------------------------------------------------------------

//{--------------------------------------------- Find 2 exons ----------------------------------------------------------------------------

char Load_Tail (char* Current_Tag, SARANGE & Tag)
{

	int Tail_Pos=0;
	char Not18=TRUE;

	if (MINIMIZE_GAP)
	{
		/* We try to extend the read until the SA gap is below Stop_Gap or if not specified till a unique hit is found */
		if(Extend_Left(Tag,Current_Tag,LH,LH-RQFACTOR,0))//left half two small for pairing
		{
		}
	}

	if (!Tag.Skip && Tag.End-Tag.Start > SAGAP_CUTOFF)
	{
		/* SA gap is too large to do pairing as it is..*/
		Tag.Level=RQFACTOR;
		Convert_To_REVSA(Tag,Current_Tag);
		Tag=Cache_SF[RQFACTOR];//save first RQFACTOR SA range...
		Tag.Conversion_Factor=Conversion_Factor-RQFACTOR;
		Not18=FALSE;
	}
	else
	{
		Convert_To_REVSA(Tag,Current_Tag);
		Tag.Level = RH;
		if (Tag.Skip) Tag.End=Tag.Start = (Conversion_Factor-revfmi->saValue[Tag.End/revfmi->saInterval]+Tag.Skip-1)-(Tag.Level);//convert unique hits to location..
		else if(Tag.Start==Tag.End) Tag.End=Tag.Start=(Conversion_Factor-BWTSaValue(revfmi,Tag.Start))-(Tag.Level);
		Tag.Conversion_Factor=Conversion_Factor-Tag.Level;
	}
	Tail_Hits_Pos[Tail_Pos++]=Tag;
	Tail_Hits_Pos[Tail_Pos].Start=0;
	return Not18;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Scan_Half_ExtendL
 *  Description:  Will scan right half for a match and try to extend it to a full match to Left..
 * =====================================================================================
 */

char Scan_Half_ExtendL(char* Current_Tag, ORIENTATION & O,unsigned Read_ID)
{
	int Head_Pos=0;
	char Short_Pairs=FALSE;//Was RQFactor-RQfactor lengths paired? 
	Unmapped_Record URecord;
	Out_Record Out;
	char Hit_Found=FALSE;

	SARANGE Tag=O.TagBH;

	if(!Load_Tail(Current_Tag,Tag))
	{
		Short_Pairs=TRUE;
	}
	//----------------------------------------------Scan Head----------------------------------------------------------------------------
	Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;
	if(Search_Forwards_Exact(Tag,Current_Tag,1,RQFACTOR)) //Left half is a match?
	{
		if (Tag.Skip) Tag.End=Tag.Start = (Conversion_Factor-revfmi->saValue[Tag.End/revfmi->saInterval]+Tag.Skip-1)-RQFACTOR;
		else if(Tag.Start==Tag.End) Tag.End=Tag.Start=(Conversion_Factor-BWTSaValue(revfmi,Tag.Start))-RQFACTOR;
		Tag.Conversion_Factor=CONVERSION_FACTOR;
		Head_Hits_Pos[Head_Pos++]=Tag;
		Head_Hits_Pos[Head_Pos].Start=0;
	//----------------------------------------------Scan Head----------------------------------------------------------------------------
		SAGAP_CUTOFF_H=SAGAP_CUTOFF_INDEX; SAGAP_CUTOFF_T= SAGAP_CUTOFF;
		Pairs[0].Head=0;
		Pair_Reads(Head_Hits_Pos,Tail_Hits_Pos);
		if(Pairs[0].Head)
		{
			Pack_Text(Current_Tag);
			int T=0,Loc,Loc1,TLoc;
			for(int i=0;T<MAX_HITS_ALLOWED && Pairs[i].Head;i++)
			{
				Tag.Start=Tag.End=Pairs[i].Head;Tag.Level=Pairs[i].HLevel;
				int MF=Extend_Location_Forward(Tag,Tag.Level,STRINGLENGTH,COUNT);//extend head further...
				Out.HLocation=Tag.Start;Out.HLength=Tag.Level;
				Tag.Start=Tag.End=Pairs[i].Tail;Tag.Level=Pairs[i].TLevel;
				int MB=Extend_Location_Backward(Tag,STRINGLENGTH-Tag.Level,0,COUNT);
				Out.TLocation=Tag.Start+1;Out.TLength=Tag.Level;
				Out.ID=Read_ID;
				Hit_Found=TRUE;
				if (Out.TLength+Out.HLength <STRINGLENGTH) continue;//reads do not extend..
				Loc=Location_To_Genome(Out.HLocation);
				Loc1=Location_To_Genome(Out.TLocation);
				if(Loc == Loc1)
				{
					Out.Type=SPLITREAD;
					Out.Strlen=STRINGLENGTH;
					Write_Sam(Current_Tag[STRINGLENGTH],Genome_Offsets[Loc].Genome,Out.HLocation,Out.TLocation,Out.HLength,Out.TLength,'S',Read_ID);
					FILE* F=Genome_Offsets[Loc].Out_File;// : Genome_Offsets[Loc].Out_FileM;
					fwrite(&Out,sizeof(Out_Record),1,F);
					if (WRITE_SPLITREAD) fprintf(Map_File, "%s\t%u\t%u\t%u\t%u\t%s\n",Genome_Offsets[Loc].Genome,Out.HLocation,Out.TLocation,Out.HLength,Out.TLength,Head.Tag_Copy);
				}
			}
			if(Hit_Found) return TRUE;
		}
	}
	if (!Hit_Found && !Short_Pairs)//check if one half matches..
	{
		if (Tail_Hits_Pos[0].End == Tail_Hits_Pos[0].Start)
		{
			Pack_Text(Current_Tag);
			Tag.Start=Tag.End=Tail_Hits_Pos[0].Start;Tag.Level=RH;
			Extend_Location_Backward(Tag,STRINGLENGTH-Tag.Level,0,COUNT_IN_EXT);
			if (Tag.Level >= STRINGLENGTH -18);// && Tag.Level <= STRINGLENGTH -10)
			{
				Tag.Start++;
				int Loc=Location_To_Genome(Tag.Start);
				URecord.Strlen=STRINGLENGTH;
				FILE* F=(O.Strand=='+') ? Genome_Offsets[Loc].UnmappedX : Genome_Offsets[Loc].UnmappedXM;
				URecord.Location=Tag.Start;URecord.Length=Tag.Level;
				Write_Sam(Current_Tag[STRINGLENGTH],Genome_Offsets[Loc].Genome,Tag.Start,0,Tag.Level,0,'T',Read_ID);
				fwrite(&URecord,sizeof(URecord),1,F);
				fwrite(Current_Tag,STRINGLENGTH,1,F);
				Hit_Found=TRUE;
				return TRUE;
			}
		}
		else //extend multiple hits..
		{
			Head_Hits_Pos[0]=O.TagBH;
			if (Head_Hits_Pos[0].End - Head_Hits_Pos[0].Start<MAX_ONE_SIDE_HITS)
			{
				Convert_To_REVSA(Head_Hits_Pos[0],Current_Tag);
				unsigned i=Head_Hits_Pos[0].Start;
				while(i<=Head_Hits_Pos[0].End)
				{
					Tag.Start=Tag.End=i;Tag.Level=RH;Tag.Skip=0;
					Tag.Start=Tag.End=Location(Tag,REVFMI)+STRINGLENGTH-RH;
					Pack_Text(Current_Tag);
					//Tag.Start=Tag.End=Tail_Hits_Pos[0].Start;Tag.Level=RH;
					Extend_Location_Backward(Tag,STRINGLENGTH-Tag.Level,0,COUNT_IN_EXT);
					if (Tag.Level >= STRINGLENGTH -18);// && Tag.Level <= STRINGLENGTH -10)
					{
						Tag.Start++;
						int Loc=Location_To_Genome(Tag.Start);
						URecord.Strlen=STRINGLENGTH;
						FILE* F=(O.Strand=='+') ? Genome_Offsets[Loc].UnmappedX : Genome_Offsets[Loc].UnmappedXM;
						URecord.Location=Tag.Start;URecord.Length=Tag.Level;
						Write_Sam(Current_Tag[STRINGLENGTH],Genome_Offsets[Loc].Genome,Tag.Start,0,Tag.Level,0,'T',Read_ID);
						fwrite(&URecord,sizeof(URecord),1,F);
						fwrite(Current_Tag,STRINGLENGTH,1,F);
						Hit_Found=TRUE;
						//return TRUE;
					}
					i++;
				}
			}
			return Hit_Found;
		}
	}

	return FALSE;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Head
 *  Description:  Loads first half of read...
 *  		  returns true if full half loaded, or returns false if quarterlength SA range loaded...
 * =====================================================================================
 */

char Load_Head (char* Current_Tag, ORIENTATION & O,SARANGE & Tag)
{

	int Head_Pos=0;
	char Not18=TRUE;
	SARANGE RQTag=O.TagFR;//save first RQFACTOR SA range...

	if (MINIMIZE_GAP)
	{
		/* We try to extend the read until the SA gap is below Stop_Gap or if not specified till a unique hit is found */
		if(Extend_Right(Tag,Current_Tag,LH+1,RH-RQFACTOR,0))//left half two small for pairing
		{
		}
	}

	if (!Tag.Skip && Tag.End-Tag.Start > SAGAP_CUTOFF)
	{
		/* SA gap is too large to do pairing as it is..*/
		Tag=RQTag;Tag.Level=RQFACTOR;
		Tag.Conversion_Factor=Conversion_Factor-Tag.Level;
		Not18=FALSE;
	}
	else
	{
		Tag.Level = LH;
		if (Tag.Skip) Tag.End=Tag.Start = (Conversion_Factor-revfmi->saValue[Tag.End/revfmi->saInterval]+Tag.Skip-1)-(Tag.Level);//convert unique hits to location..
		else if(Tag.Start==Tag.End) Tag.End=Tag.Start=(Conversion_Factor-BWTSaValue(revfmi,Tag.Start))-(Tag.Level);
		Tag.Conversion_Factor=Conversion_Factor-Tag.Level;
	}
	Head_Hits_Pos[Head_Pos++]=Tag;
	Head_Hits_Pos[Head_Pos].Start=0;
	return Not18;
}

char Load_HeadX (char* Current_Tag, ORIENTATION & O,SARANGE & Tag)
{

	/* SA gap is too large to do pairing as it is..*/
	Tag=O.TagFR;Tag.Level=RQFACTOR;
	Tag.Conversion_Factor=Conversion_Factor-Tag.Level;
	if (!Tag.Start) return FALSE;

	if (Tag.Skip) Tag.End=Tag.Start = (Conversion_Factor-revfmi->saValue[Tag.End/revfmi->saInterval]+Tag.Skip-1)-(Tag.Level);//convert unique hits to location..
	else if(Tag.Start==Tag.End) Tag.End=Tag.Start=(Conversion_Factor-BWTSaValue(revfmi,Tag.Start))-(Tag.Level);
	Tag.Conversion_Factor=Conversion_Factor-Tag.Level;

	Head_Hits_Pos[0]=Tag;
	Head_Hits_Pos[1].Start=0;
	return TRUE;
}

char Scan_Half_ExtendR(char* Current_Tag, ORIENTATION & O,unsigned Read_ID)
{
	int Tail_Pos=0;
	SARANGE Tag=O.TagFH;
	char Short_Pairs=FALSE;//Was RQFactor-RQfactor lengths paired? 
	Out_Record Out;
	Unmapped_Record URecord;
	int Hit_Found=FALSE;

	
	if(!Load_Head(Current_Tag,O,Tag))
	{
		Short_Pairs=TRUE;
	}
	//----------------------------------------------Scan Tail----------------------------------------------------------------------------
	Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;
	if(Search_Forwards_Exact(Tag,Current_Tag,STRINGLENGTH-RQFACTOR,RQFACTOR)) //Scan for last RQFACTOR...
	{
		if (Tag.Skip) Tag.End=Tag.Start = (Conversion_Factor-revfmi->saValue[Tag.End/revfmi->saInterval]+Tag.Skip-1)-RQFACTOR;
		else if(Tag.Start==Tag.End) Tag.End=Tag.Start=(Conversion_Factor-BWTSaValue(revfmi,Tag.Start))-RQFACTOR;
		Tag.Conversion_Factor=CONVERSION_FACTOR;
		Tail_Hits_Pos[Tail_Pos++]=Tag;
		Tail_Hits_Pos[Tail_Pos].Start=0;
	//----------------------------------------------Scan Tail----------------------------------------------------------------------------
		SAGAP_CUTOFF_T=SAGAP_CUTOFF_INDEX; SAGAP_CUTOFF_H= SAGAP_CUTOFF;
		Pairs[0].Head=0;
		Pair_Reads(Head_Hits_Pos,Tail_Hits_Pos);
		if(Pairs[0].Head)
		{
			Pack_Text(Current_Tag);
			int Loc,Loc1,TLoc,T=0;
			for(int i=0;T<MAX_HITS_ALLOWED && Pairs[i].Head;i++)
			{
				Tag.Start=Tag.End=Pairs[i].Head;Tag.Level=Pairs[i].HLevel;
				Extend_Location_Forward(Tag,Tag.Level,STRINGLENGTH,COUNT);//extend head further...
				if (Short_Pairs && Tag.Level < LH) continue;
				Out.HLocation=Tag.Start;Out.HLength=Tag.Level;
				Tag.Start=Tag.End=Pairs[i].Tail+1;Tag.Level=Pairs[i].TLevel;
				Extend_Location_Backward(Tag,STRINGLENGTH-Tag.Level,0,COUNT);
				Out.TLocation=Tag.Start+1;Out.TLength=Tag.Level;
				if (Out.TLength+Out.HLength <STRINGLENGTH) continue;//reads do not extend..
				Hit_Found=TRUE;
				Out.Type=UNEXTENDED;
				Out.ID=Read_ID;
				Loc=Location_To_Genome(Out.HLocation);
				Loc1=Location_To_Genome(Out.TLocation);
				if (Loc == Loc1)//in the same chromosome?
				{
					Write_Sam(Current_Tag[STRINGLENGTH],Genome_Offsets[Loc].Genome,Out.HLocation,Out.TLocation,Out.HLength,Out.TLength,'S',Read_ID);
					Out.Type=SPLITREAD;
					Out.Strlen=STRINGLENGTH;
					FILE* F= Genome_Offsets[Loc].Out_File;// : Genome_Offsets[Loc].Out_FileM;
					fwrite(&Out,sizeof(Out),1,F);
					if (WRITE_SPLITREAD) fprintf(Map_File, "%s\t%u\t%u\t%u\t%u\t%s\n",Genome_Offsets[Loc].Genome,Out.HLocation,Out.TLocation,Out.HLength,Out.TLength,Head.Tag_Copy);
				}
			}
			if(Hit_Found) return TRUE;
		}
	}
	if (!Hit_Found && !Short_Pairs)//check if one half matches..
	{
		if (Head_Hits_Pos[0].End == Head_Hits_Pos[0].Start)//extend single hit
		{
			Pack_Text(Current_Tag);
			Tag.Start=Tag.End=Head_Hits_Pos[0].Start;Tag.Level=LH;
			Extend_Location_Forward(Tag,Tag.Level,STRINGLENGTH,COUNT);//extend head further...
			if (Tag.Level >= STRINGLENGTH -18)// && Tag.Level <= STRINGLENGTH -10)
			{
				int Loc=Location_To_Genome(Tag.Start);
				FILE* F=(O.Strand=='+') ? Genome_Offsets[Loc].Unmapped : Genome_Offsets[Loc].UnmappedM;
				URecord.Location=Tag.Start;URecord.Length=Tag.Level;
				URecord.Strlen=STRINGLENGTH;
				Write_Sam(Current_Tag[STRINGLENGTH],Genome_Offsets[Loc].Genome,Tag.Start,0,Tag.Level,0,'T',Read_ID);
				fwrite(&URecord,sizeof(URecord),1,F);
				fwrite(Current_Tag,STRINGLENGTH,1,F);
				return (Hit_Found=TRUE);
			}
		}
		else //extend multiple hits..
		{
			Head_Hits_Pos[0]=O.TagFH;
			if (Head_Hits_Pos[0].End - Head_Hits_Pos[0].Start<MAX_ONE_SIDE_HITS)
			{
				unsigned i=Head_Hits_Pos[0].Start;
				while(i<=Head_Hits_Pos[0].End)
				{
					Tag.Start=Tag.End=i;Tag.Level=LH;Tag.Skip=0;
					Tag.Start=Tag.End=Location(Tag,REVFMI)+STRINGLENGTH-LH;
					Pack_Text(Current_Tag);
					Extend_Location_Forward(Tag,Tag.Level,STRINGLENGTH,COUNT);//extend head further...
					if (Tag.Level >= STRINGLENGTH -18 );//&& Tag.Level <= STRINGLENGTH -10)
					{
						int Loc=Location_To_Genome(Tag.Start);
						FILE* F=(O.Strand=='+') ? Genome_Offsets[Loc].Unmapped : Genome_Offsets[Loc].UnmappedM;
						URecord.Location=Tag.Start;URecord.Length=Tag.Level;
						URecord.Strlen=STRINGLENGTH;
						Write_Sam(Current_Tag[STRINGLENGTH],Genome_Offsets[Loc].Genome,Tag.Start,0,Tag.Level,0,'T',Read_ID);
						fwrite(&URecord,sizeof(URecord),1,F);
						fwrite(Current_Tag,STRINGLENGTH,1,F);
					}
					i++;
				}
			}
		}
	}
	return FALSE;
}


//}--------------------------------------------- Find 2 exons ----------------------------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Pack_Text
 *  Description:  Binary encode read...
 * =====================================================================================
 */
void Pack_Text(char* Current_Tag)
{
	int i=0,j=0,k=0;
	char Temp_Char=0;
	while(i< STRINGLENGTH)
	{
		Temp_Char <<= 2;
		Temp_Char |= Current_Tag[i];
		if (++j == 4)//one char packed...
		{
			Packed[k]=Temp_Char;
			Temp_Char=0;j=0;
			k++;
		}
		i++;
	}
	if(k<((STRINGLENGTH-1)/4)) Packed[k]=Temp_Char;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Extend_Location
 *  Description:  Search and extend hit in the original location...
 *  		  Starts search from Start and until tag [STRINGLENGTH-Start .. STRINGLENGTH-Start-Stringlength] has been searched...
 *  		  allow Count mismatches.. Tag.Start=Tag.End=location of search...
 * =====================================================================================
 */
int Extend_Location_Backward(SARANGE & Tag,int Start,int End,int Count)
{
	unsigned Source,Destination,Xor;
	unsigned Org,Piece_Length;
	int Shift,Mismatches=0,StringLength;
	char Now;
	int Mis_Level=0;

	Org=Tag.Start;
	StringLength=Start-End;
	Start++;
	Tag.Level=1;

	for(;;)
	{

		Source=(BSWAP(*(UINT*)(Original_Text+(Org/4)-BYTES_IN_UNIT+1)) >> (8-(Org % 4) * 2));//>>30;;
		Destination=(BSWAP(*(UINT*)(Packed-BYTES_IN_UNIT+1+(Start-Tag.Level)/4)) >> (8-((Start-Tag.Level) % 4) * 2));//>>30;
		Piece_Length=StringLength-Tag.Level;

		if ( Piece_Length > MULTISTRINGLENTH ) Piece_Length=MULTISTRINGLENTH;
		Shift=INTLENGTH-2*(Piece_Length);
		Source <<= Shift;Destination <<= Shift;
		Source >>= Shift;Destination >>= Shift;
		Source &= SEARCH_MASKB;
		Destination &= SEARCH_MASKB;

		if (!(Xor =Destination ^ Source))
		{	
			Org -= Piece_Length;
			Tag.Level += Piece_Length; 
			//if(Tag.Level> StringLength)
			if(Tag.Level == StringLength)
			{
				Tag.Start=Tag.End=Org;
				Tag.Level+=STRINGLENGTH-Start-1;
				return Tag.Level;	
			}
			else {continue;}
		} 
		else//mismatch...
		{
			//check if only one mismatch...
			int First_One = BSR(Xor);
			First_One = First_One /2 ;
			Org -= First_One+1;
			Tag.Level += First_One;
			Now = Source >> 2*First_One;
			Now &= 3;

			//Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			//Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level-1);
			//Tag.Mismatches++;
			Mismatches++;if(!Mis_Level) Mis_Level=Tag.Level+STRINGLENGTH-Start;
			Tag.Level++;

			if (Mismatches<=Count) 
			{
				if(Tag.Level == StringLength)
				{
					Tag.Start=Tag.End=Org;
					Tag.Level+=STRINGLENGTH-Start-1;
					return Mis_Level;	
				}
				else {continue;}
			} 
			else//log 2 mismatches 
			{
				Tag.Start=Tag.End=Org;
				Tag.Level+=STRINGLENGTH-Start-1;
				return Mis_Level;	
				//return Mis_Level;
			}
		}

	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Extend_Location
 *  Description:  Search and extend hit in the original location...
 *  		  Starts search from Start and until tag [1..StringLength] has been searched...
 *  		  allow Count mismatches.. Tag.Start=Tag.End=location of search...
 * =====================================================================================
 */
int Extend_Location_Forward(SARANGE & Tag,int Start,int StringLength,int Count)
{
	unsigned Source,Destination,Xor;
	unsigned Org,Piece_Length;
	int Shift,Mismatches=0;
	char Now;
	int Mis_Level=0;

	Org=Tag.Start+Start;

	for(;;)
	{

		Source=(BSWAP(*(UINT*)(Original_Text+Org/4)) << ((Org % 4) * 2));//>>30;;
		Destination=(BSWAP(*(UINT*)(Packed+(Start)/4)) << (((Start) % 4) * 2));//>>30;
		Piece_Length=StringLength-Start+1;

		if ( Piece_Length > MULTISTRINGLENTH ) Piece_Length=MULTISTRINGLENTH;
		Shift=INTLENGTH-2*(Piece_Length);
		Source >>= Shift;Destination >>= Shift;
		Source <<= Shift;Destination <<= Shift;
		Source &= SEARCH_MASK;
		Destination &= SEARCH_MASK;

		if (!(Xor =Destination ^ Source))
		{
			Org += Piece_Length-1;
			Start += Piece_Length-1;
			if(Start== StringLength)
			{
				//Tag.Start=Tag.End=Org;
				Tag.Level=Start;
				return Tag.Level;	
			}
			else {Org++;Start++;continue;}
		} 
		else//mismatch...
		{
			//check if only one mismatch...
			int First_One = BSF(Xor);
			First_One = First_One /2 ;
			Org += First_One;
			Start += First_One;
			Now = Source << 2*First_One;
			Now >>= SEARCH_ISOLATE_SHIFT;

			//Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			//Tag.Mismatch_Pos[Tag.Mismatches]=(Start);
			//Tag.Mismatches++;
			Mismatches++;if(!Mis_Level) Mis_Level=Start;
			if (Mismatches<=Count) 
			{
				if(Start == StringLength)
				{
					//Tag.Start=Tag.End=Org;
					Tag.Level=Start;
					return Mis_Level;	
				}
				else {Org++;Start++;continue;}
			} 
			else
			{
				//Tag.Start=Tag.End=Org;
				Tag.Level=Start;
				return Mis_Level;
			}

		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Init 
 *  Description:  Load indexes, Files and initialise variables.. 
 * =====================================================================================
 */

void Init()
{


	Open_Files();
	if (MAPMODE)
	{
		Load_Range_Index(INDFILE,BLKFILE,Range_Index);
		SA_Index=Range_Index.SA_Index;
		SA_Blocks=Range_Index.SA_Blocks;
		COMPRESS=Range_Index.COMPRESS;
		Hits=Range_Index.Hits;
	}


	Original_File=File_Open(BINFILE,"rb");
	if(MAPMODE) SAM=File_Open("Mapped.sam","w");
	Original_Text=(unsigned char*) malloc(Get_File_Size(Original_File));
	fread(Original_Text,Get_File_Size(Original_File),1,Original_File);

	Detect_Input(TAG_COPY_LEN, STRINGLENGTH, NORMAL_TAGS, PAIRING_TYPE, FILETYPE,PAIR_LENGTH_RIGHT, PAIR_LENGTH_LEFT);
	if (NORMAL_TAGS) printf("Single end reads...\n"); else printf("Paired end reads...\n");
	if(MAPMODE)
	{
		if(MINIMUMINDEX) fwfmi= Load_Indexes(BWTFILE,OCCFILE,NULL); else fwfmi= Load_Indexes(BWTFILE,OCCFILE,SAFILE);
		revfmi= Load_Indexes(REVBWTINDEX,REVOCCFILE,REVSAFILE);
		fwfmi->saInterval=revfmi->saInterval;

		SOURCELENGTH = revfmi->textLength;
		Conversion_Factor=revfmi->textLength;//+1;
		CONVERSION_FACTOR=revfmi->textLength-RQFACTOR;//+1;
	}

	/*if(!NORMAL_TAGS)
	{
		PRLH=PAIR_LENGTH_RIGHT/2;//calculate tag portions...
		if ((PAIR_LENGTH_RIGHT % 2)) {PRLH++;}	
		PRRH=PAIR_LENGTH_RIGHT-PRLH;
		PLLH=PAIR_LENGTH_LEFT/2;//calculate tag portions...
		if ((PAIR_LENGTH_LEFT % 2)) {PLLH++;}	
		PLRH=PAIR_LENGTH_LEFT-PLLH;
	}
	else*/
	{
		STRINGLENGTHO=STRINGLENGTH;
		LH=STRINGLENGTH/2;//calculate tag portions...
		if ((STRINGLENGTH % 2)) {LH++;}	
		RH=STRINGLENGTH-LH;
		RHQL=RH/2;RHQR=RH-RHQL;
	}
	Char_To_Code['N']=0;Char_To_Code['n']=0;Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt

	Allocate_Memory();

	BigBed=File_Open(WIGGLEFILE,"w");
	Junctions_File=File_Open(JUNCTIONFILE,"w");
	fprintf(BigBed,"track type=bedGraph name=\"Rawbin Coverage\"\n");
//strt with 256
	Int_To_Junc[0]= "+";
	Int_To_Junc[1]= "GT-AG";
	Int_To_Junc[2]= "GC-AG";
	Int_To_Junc[3]= "AC-AT";

	Int_To_Junc[4]= "CT-AC";
	Int_To_Junc[5]= "CT-GC";
	Int_To_Junc[6]= "AT-GT";
	Int_To_Junc[7]= "-";

}


void Allocate_Memory()
{
	Buffer=(Out_Record*)malloc(BUFSIZE*sizeof(Out_Record));

	int STRINGLENGTH =36;
	int Max_Allocate=1;
	int Max_Limit=5;

	Head_Hits_Neg=(SARANGE*)malloc((MAX_HITS_TO_STORE+1)*sizeof(SARANGE));
	Head_Hits_Pos=(SARANGE*)malloc((MAX_HITS_TO_STORE+1)*sizeof(SARANGE));
	Tail_Hits_Neg=(SARANGE*)malloc((MAX_HITS_TO_STORE+1)*sizeof(SARANGE));
	Tail_Hits_Pos=(SARANGE*)malloc((MAX_HITS_TO_STORE+1)*sizeof(SARANGE));
	if (!(Pairs=(PAIR*)malloc(sizeof(PAIR)*(MAX_HITS_TO_STORE+10)))) {printf("Allocate_Memory():malloc error...\n");exit(0);}

	FSSStack=(SARANGE*)malloc(sizeof(SARANGE)*2*MAXSTRINGLEN);	
	BMStack=(SARANGE*)malloc(sizeof(SARANGE)*2*MAXSTRINGLEN);	
	BMHStack=(SARANGE*)malloc(sizeof(SARANGE)*2*MAXSTRINGLEN);	

	for (int i=0;i<Max_Limit-1;i++) Max_Allocate=Max_Allocate*STRINGLENGTH;

	ARRAY_BOUND =sizeof(SARANGE)*4*Max_Allocate;
        END_BOUND=sizeof(SARANGE)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH;
	Mismatches_ForwardP=(SARANGE*)malloc(ARRAY_BOUND);
	Mismatches_ForwardC=(SARANGE*)malloc(ARRAY_BOUND);
	Mismatches_BackwardP=(SARANGE*)malloc(ARRAY_BOUND);
	Mismatches_BackwardC=(SARANGE*)malloc(ARRAY_BOUND);
	Two_Mismatches_At_EndP=(SARANGE*)malloc(END_BOUND);	
	Two_Mismatches_At_EndC=(SARANGE*)malloc(END_BOUND);	
	Two_Mismatches_At_End_ForwardP=(SARANGE*)malloc(END_BOUND);	
	Two_Mismatches_At_End_ForwardC=(SARANGE*)malloc(END_BOUND);	

	Can_Acceptor=(State*)malloc(sizeof(State)*OVERHANG);
	Can_Donor=(State*)malloc(sizeof(State)*OVERHANG);

}

//{-----------------------------  BACKWARD SEARCH ROUTINE  -------------------------------------------------/

void Convert_To_REVSA(SARANGE & Tag, char* Current_Tag)
{
			//int StringLength=Tag.Level+RH-1;
			int StringLength=Tag.Level;
			Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;
#ifdef DEBUG
			if(Search_Forwards_Exact(Tag,Current_Tag,STRINGLENGTH-StringLength+1,StringLength))//convert to reverse SA range...
				printf("Convert_To_REVSA()!: Invalid SA range...\n");
#else			
			Search_Forwards_Exact(Tag,Current_Tag,STRINGLENGTH-StringLength+1,StringLength);//convert to reverse SA range...
#endif
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Extend_Right
 *  Description:  Extend the SA range in tag 1, until fail ,string length amount scanned
 *  		  or SAGap < Stop_Gap.
 * =====================================================================================
 */
char Extend_Right(SARANGE & Tag,char* Current_Tag,int Start,int StringLength,int Stop_Gap)
{
	unsigned Index,First,Last;
	SARANGE Temp;//temporary tag to save last tag details before failure...

	char Now;
	Tag.Level=1;
	Start -= 2;//accomodate 1 based offset.. Start is 1 based and so is tag Level
	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}

			for(;;)
			{
				Index=Tag.Start;
				if (Index >= revfmi->inverseSa0) Index--;//adjust for missing $
				Now=revfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
				if (Current_Tag[Start+Tag.Level] == Now)
				{
					Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;
					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}

					if(Tag.Level== StringLength)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						return TRUE; 	
					}
					else {Tag.Level++;continue;}
				} 
				else//mismatch...
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					return FALSE;	
				}
			}
		}
		else//SA range has sevaral possible hits... 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= revfmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned Pos=First;Pos<=Last;Pos++)
				{
					Now=revfmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Tag.Level+Start];
				if (Branch_Characters[Now])//we have a match... 
				{
					Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				}
				else//mismatch..
				{
					return FALSE;
				}
			} 
			else
			{
				Temp=Tag;
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,revfmi);
				if (!Tag.Start) {Tag=Temp;return FALSE;}
			}

			if(Tag.Level== StringLength)
			{
				return TRUE;
			}
			else 
			{
				Tag.Level++;
				if(Stop_Gap && (Tag.End-Tag.Start) < Stop_Gap) return FALSE;
				continue;
			}
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Extend_Left
 *  Description:  Extend the SA range in tag 1, until fail ,string length amount scanned
 *  		  or SAGap < Stop_Gap.
 *  		  returns TRUE if full stringlength scanned..
 * =====================================================================================
 */
char Extend_Left(SARANGE & Tag,char* Current_Tag,int Start,int StringLength,int Stop_Gap)
{

	unsigned Index,First,Last;
	char Now;
	SARANGE Temp;
	Tag.Level=1;

	if (!Tag.Start) return FALSE;
	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;
				Tag.End=Tag.Start;
			}
			for(;;)
			{
				Index=Tag.Start;
				if (Index >= fwfmi->inverseSa0) Index--;//adjust for missing $
				Now=fwfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
				if (Current_Tag[Start-Tag.Level] == Now)
				{
					Tag.Start = fwfmi->cumulativeFreq[Now] + BWTOccValue(fwfmi, Tag.Start, Now) + 1;

					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}

					if(Tag.Level== StringLength)//no need to print as we are still halfway..
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						return TRUE;	
					}
					else {Tag.Level++;continue;}
				} 
				else
				{
					return FALSE;	
				}
			}
		}
		else 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fwfmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fwfmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Start-Tag.Level];
				if (Branch_Characters[Now])//we have a match... 
				{

					Tag.Start = fwfmi->cumulativeFreq[Now] + BWTOccValue(fwfmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				
				}
				else return FALSE;
			}
			else 
			{
				Temp=Tag;
				Get_SARange_Fast(Current_Tag[Start-Tag.Level],Tag,fwfmi);
				if (!Tag.Start) {Tag=Temp;return FALSE;}
			}

			if(Tag.Level== StringLength) return TRUE;
			else 
			{
				Tag.Level++;
				if(Stop_Gap && (Tag.End-Tag.Start)< Stop_Gap) return FALSE;
				continue;
			}

		}
	}
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Backwards_Exact
 *  Description:  backward seach for exact occurence for string in Current_Tag in fwfmi
 *  		  Start is the 1-indexed location of string position in Current_Tag
 *  		  String length is the length of substring to search...
 * =====================================================================================
 */

char Search_Backwards_Exact(SARANGE & Tag,char* Current_Tag,int Start,int StringLength)
{

	unsigned Index,First,Last;
	char Now;
	SARANGE Temp;
	Tag.Level=1;

	if (!Tag.Start) return FALSE;
	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;
				Tag.End=Tag.Start;
			}
			for(;;)
			{
				Index=Tag.Start;
				if (Index >= fwfmi->inverseSa0) Index--;//adjust for missing $
				Now=fwfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
				if (Current_Tag[Start-Tag.Level] == Now)
				{
					Tag.Start = fwfmi->cumulativeFreq[Now] + BWTOccValue(fwfmi, Tag.Start, Now) + 1;

					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}
					else {Tag.End=Tag.Start;}

					Cache_SB[Tag.Level]=Tag;
					if(Tag.Level== StringLength)//no need to print as we are still halfway..
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						return TRUE;	
					}
					else {Tag.Level++;continue;}
				} 
				else
				{
					Tag.Level--;
					return FALSE;	
				}
			}
		}
		else 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fwfmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fwfmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Start-Tag.Level];
				if (Branch_Characters[Now])//we have a match... 
				{

					Tag.Start = fwfmi->cumulativeFreq[Now] + BWTOccValue(fwfmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				
				}
				else {Tag.Level--;return FALSE;}
			}
			else 
			{
				Temp=Tag;
				Get_SARange_Fast(Current_Tag[Start-Tag.Level],Tag,fwfmi);
				if (!Tag.Start) {Tag=Temp;return FALSE;}
			}

			Cache_SB[Tag.Level]=Tag;
			if(Tag.Level== StringLength) return TRUE;
			else {Tag.Level++;continue;}

		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Forwards_Exact
 *  Description:  forward seach for exact occurence for string in Current_Tag in revfmi
 *  		  Start is the 1-indexed location of string position in Current_Tag
 *  		  String length is the length of substring to search...
 *  		  Start+Tag.Length-2 last position...
 *  		  return TRUE if full StringLength scanned.
 * =====================================================================================
 */

char Search_Forwards_Exact(SARANGE & Tag,char* Current_Tag, int Start,int StringLength)
{
	unsigned Index,First,Last;
	SARANGE Temp;//temporary tag to save last tag details before failure...

	char Now;
	Tag.Level=1;
	Start -= 2;//accomodate 1 based offset.. Start is 1 based and so is tag Level
	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}

			for(;;)
			{
				Index=Tag.Start;
				if (Index >= revfmi->inverseSa0) Index--;//adjust for missing $
				Now=revfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
				if (Current_Tag[Start+Tag.Level] == Now)
				{
					Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;
					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}
					else {Tag.End=Tag.Start;}

					Cache_SF[Tag.Level]=Tag;
					if(Tag.Level== StringLength)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						return TRUE; 	
					}
					else {Tag.Level++;continue;}
				} 
				else//mismatch...
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level--;
					return FALSE;	
				}
			}
		}
		else//SA range has sevaral possible hits... 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= revfmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned Pos=First;Pos<=Last;Pos++)
				{
					Now=revfmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Tag.Level+Start];
				if (Branch_Characters[Now])//we have a match... 
				{
					Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				}
				else//mismatch..
				{
					return FALSE;
				}
			} 
			else
			{
				Temp=Tag;
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,revfmi);
				if (!Tag.Start) {Tag=Temp;return FALSE;}
			}

			Cache_SF[Tag.Level]=Tag;
			if(Tag.Level== StringLength)
			{
				return TRUE;
			}
			else {Tag.Level++;continue;}

		}
	}
}
//}-----------------------------  BACKWARD SEARCH ROUTINE  -------------------------------------------------/

//{----------------------------------- FILE HANDLING ---------------------------------------------------------


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_Tag 
 *  Description:  Read next tag and store details in Head.. 
 * =====================================================================================
 */

char Read_Tag(READ & Head,READ & Tail)
{
	char * Current_Tag;
	static int Random_Pointer=0;

	//if(NORMAL_TAGS)
	//{
	Current_Tag=Head.Tag;
	if (gzgets(Input_File,Head.Description,MAXDES)!=0)// read a tag...
	{
		gzgets(Input_File,Head.Tag,MAXDES);//tag
		if (FILETYPE == FQ)
		{
			gzgets(Input_File,Head.Plus,MAXTAG);//plus
			gzgets(Input_File,Head.Quality,MAXTAG);//phred
		}
		strcpy(Head.Tag_Copy,Head.Tag);
		Head.NCount=0;int j=0;
		for (unsigned i=0;i<=STRINGLENGTH-1;i++)
		{
			if (Current_Tag[i] == 'n' || Current_Tag[i]=='N')
			{
				Head.N[j++]=i;Head.NLocations[i]=TRUE;Head.NCount++;
				Current_Tag[i]=Random_Array[Random_Pointer++];Head.N[j++]=Current_Tag[i];
				if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
			}
			else Head.NLocations[i]=FALSE;
			Current_Tag[i]=Char_To_Code[Current_Tag[i]];
			Head.Complement[STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];
		}
		Current_Tag[STRINGLENGTH]='+';
		Head.Complement[STRINGLENGTH]='-';
	} else return FALSE;
	if(NORMAL_TAGS) return TRUE;

	Current_Tag=Tail.Tag;
	if (gzgets(Mate_File,Tail.Description,MAXDES)!=0)// read a tag...
	{
		gzgets(Mate_File,Tail.Tag,MAXDES);//tag
		if (FILETYPE == FQ)
		{
			gzgets(Mate_File,Tail.Plus,MAXTAG);//plus
			gzgets(Mate_File,Tail.Quality,MAXTAG);//phred
		}
		strcpy(Tail.Tag_Copy,Tail.Tag);
		Tail.NCount=0;int j=0;
		for (unsigned i=0;i<=STRINGLENGTH-1;i++)
		{
			if (Current_Tag[i] == 'n' || Current_Tag[i]=='N')
			{
				Tail.N[j++]=i;Tail.NLocations[i]=TRUE;Tail.NCount++;
				Current_Tag[i]=Random_Array[Random_Pointer++];Tail.N[j++]=Current_Tag[i];
				if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
			}
			else Tail.NLocations[i]=FALSE;
			Current_Tag[i]=Char_To_Code[Current_Tag[i]];
			Tail.Complement[STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];
		}
		Current_Tag[STRINGLENGTH]='+';
		Tail.Complement[STRINGLENGTH]='-';
		Current_Tag=Head.Tag;
		return TRUE;
	} else {printf("Read_Tag():Unpaired read...!\n");exit(0);};
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Detect_Input
 *  Description:  Detects file type and lengths of input file(s)
 * =====================================================================================
 */
void Detect_Input(int & TAG_COPY_LEN,int & STRINGLENGTH,char & NORMAL_TAGS, char & PAIRING_TYPE,char & FILETYPE,int & PAIR_LENGTH_RIGHT,int & PAIR_LENGTH_LEFT)
{
	char Description[MAXDES+1];
	char Current_Tag[MAXTAG+1];
	char Quality[MAXTAG+1];

	if (gzgets(Input_File,Description,MAXDES)!=0)//Measure tag length
	{
		gzgets(Input_File,Current_Tag,MAXTAG);
		for(TAG_COPY_LEN=0;Current_Tag[TAG_COPY_LEN]!='\n' && Current_Tag[TAG_COPY_LEN]!='\r' && Current_Tag[TAG_COPY_LEN]!=0;TAG_COPY_LEN++);//TAG_COPY_LEN++;
		if(Patternfile_Count)
		{
			gzgets(Mate_File,Description,MAXDES);
			Current_Tag[TAG_COPY_LEN++]='\t';gzgets(Mate_File,Current_Tag+TAG_COPY_LEN,MAXTAG);
			for(TAG_COPY_LEN=0;Current_Tag[TAG_COPY_LEN]!='\n' && Current_Tag[TAG_COPY_LEN]!='\r' && Current_Tag[TAG_COPY_LEN]!=0;TAG_COPY_LEN++);//TAG_COPY_LEN++;
		}
		for(STRINGLENGTH=0;Current_Tag[STRINGLENGTH]!='\n' && Current_Tag[STRINGLENGTH]!='\r' && Current_Tag[STRINGLENGTH]!=0 && Current_Tag[STRINGLENGTH]!=PAIR_END_SEPERATOR;STRINGLENGTH++);
		if(Current_Tag[STRINGLENGTH]==PAIR_END_SEPERATOR) 
		{
			NORMAL_TAGS=FALSE;//we have pair ended tags..
			if(Patternfile_Count) {PAIRING_TYPE=TWOFILE;}else {PAIRING_TYPE=TAB;}
			PAIR_LENGTH_LEFT=STRINGLENGTH;

			for(PAIR_LENGTH_RIGHT=0;Current_Tag[STRINGLENGTH+1+PAIR_LENGTH_RIGHT]!='\n' && Current_Tag[STRINGLENGTH+1+PAIR_LENGTH_RIGHT]!='\r' && Current_Tag[STRINGLENGTH+1+PAIR_LENGTH_RIGHT]!=0;PAIR_LENGTH_RIGHT++);
			gzgets(Input_File,Quality,MAXTAG);//plus
			if (Quality[0]=='>') FILETYPE=FA;else FILETYPE=FQ;
			if (FILETYPE == FQ && Quality[0] != '+' && Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
		}
		else
		{
			NORMAL_TAGS=TRUE;
			gzgets(Input_File,Quality,MAXTAG);//plus
			if (Quality[0]=='>') FILETYPE=FA;else FILETYPE=FQ;
			if (FILETYPE == FQ && Quality[0] != '+' && Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
			gzgets(Input_File,Quality,MAXTAG);//phred
		}

		fseek(Input_FileO, 0L, SEEK_END);
		File_Size = ftello64(Input_FileO);

		gzseek(Input_File,0,SEEK_SET);//go top
		gzseek(Mate_File,0,SEEK_SET);//go top
	}
}

void Open_Files()
{
	Input_File=File_OpenZ(PATTERNFILE,"r");//Load tags
	if(Patternfile_Count) Mate_File=File_OpenZ(PATTERNFILE1,"r");//Load tags
	if (MAPMODE && WRITE_SPLITREAD) Map_File=File_Open(MAPFILE,"w");
	gz_stream *s=(gz_stream*)Input_File;
	Input_FileO=s->file;
	Junc_Log=File_Open("Junc.ID","w");
}

unsigned Get_File_Size(FILE* File)
{
	fseek (File , 0 , SEEK_END);
	unsigned Size = ftell (File);
	rewind (File);
	return Size;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  File_Open
 *  Description:  Open a file:
 *  Mode - "w" create/Write text from top    - "wb" Write/Binary  -"w+" open text read/write            -"w+b" same as prev/Binary
 *         "r" Read text from top            - "rb" Read/Binary   -"r+" open text read/write/nocreate   -"r+b" same as prev/binary
 *       - "a" text append/write                                  -"a+" text open file read/append      -"a+b" open binary read/append
 *
 * =====================================================================================
 */
FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File_Open:File %s Cannot be opened! \n",File_Name);
		exit(1);
	}
	else return Handle;
}

gzFile File_OpenZ(const char* File_Name,const char* Mode)
{
	gzFile Handle;
	Handle=gzopen(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File_OpenZ:File %s Cannot be opened ....\n",File_Name);
		exit(1);
	}
	return Handle;
}
//}----------------------------------- FILE HANDLING ---------------------------------------------------------

//{----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Location
 *  Description:  return the raw coordinate in the index. for revfmi, need to adjust with 
 *  		  the length of the scanned tag...
 * =====================================================================================
 */
unsigned Location(SARANGE & Tag, char Index)
{

	//Location(Tag,REVFMI)-Tag.Level to get correct location..
	if(Index)//in revfmi...
	{
		if (Tag.Skip) return Conversion_Factor-STRINGLENGTH-revfmi->saValue[Tag.End/revfmi->saInterval]+Tag.Skip-1;
		else  return Conversion_Factor-STRINGLENGTH-BWTSaValue(revfmi,Tag.Start);
	}
	else
	{
		if (Tag.Skip) return fwfmi->saValue[Tag.End/fwfmi->saInterval]-Tag.Skip+1;
		else return BWTSaValue(fwfmi,Tag.Start);
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Print_LocationX(SARANGE & Tag,char Index)
 *  Description:  Print a hit correspoding to Tag[0] in Index to file..
 * =====================================================================================
 */
void Print_LocationX(SARANGE & Tag,char Index)
{
	unsigned L=Location(Tag,Index);
	unsigned Loc=Location_To_Genome(L);
	Out_Record Out;

	Out.HLength=Out.Strlen=STRINGLENGTH;
	Out.HLocation=L;
	Out.TLength=0;
	Out.TLocation=0;
	Out.Type=FULLREAD;
	fwrite(&Out,sizeof(Out),1,Genome_Offsets[Loc].Out_File);
	Write_Sam(Current_Tag[STRINGLENGTH],Genome_Offsets[Loc].Genome,L,0,STRINGLENGTH,0,'M',0);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Indexes 
 *  Description:  Opens FM index 
 * =====================================================================================
 */

BWT* Load_Indexes(char *BWTINDEX,char *OCCFILE, char *SAFILE)
{
        int PoolSize = 524288;
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);
	return BWTLoad(mmPool, BWTINDEX, OCCFILE, SAFILE, NULL, NULL, NULL);//Load FM index
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  UnLoad_Indexes 
 *  Description:  Opens FM index 
 * =====================================================================================
 */

void UnLoad_Indexes()
{
	BWTFree(mmPool,fwfmi);
	BWTFree(mmPool,revfmi);
	MMPoolFree(mmPool);
	free(Range_Index.SA_Index);
	free(Range_Index.SA_Blocks);
}

void Get_SARange_Fast( char New_Char, SARANGE & Range,BWT *fmi)
{
	Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.Start, New_Char) + 1;
	Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.End+1, New_Char);
	if (Range.End<Range.Start) 
	{
		Range.Start=0;
	}
}
//}----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------

//{-----------------------------  Parse Command Line  -------------------------------------------------
void Parse_Command_line(int argc, char* argv[])
{
	int Current_Option=0;
	char Short_Options[] ="jM:bhq:t:g:G:n:N:o:w:mprO::";//allowed options....
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
" --query | -n <number>\t\t number of mismatches in non spliced junctions...\n"
" --query | -N <number>\t\t number of mismatches in splices...\n"
" --query | -b \t\t build files from refgene\n"
" --query | -j \t\t output all junctions...\n"
;

	if(argc == 1) {printf("%s \n",Help_String);exit(0);}
	Source=(char*)malloc(sizeof(char)*6000);//create space for file names...
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
			case 'j':
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
			case 't':
				MAX_TAGS_TO_PROCESS=atoi(optarg);
				break;
			case 'G':
				EXONGAP=atoi(optarg);
				break;
			case 'n':
				MIS_IN_INITMAP=atoi(optarg);
				break;
			case 'N':
				COUNT=atoi(optarg);
				break;
			case 'r':
				USEREFGENE=TRUE;
				break;
			case 'q':

				if(!Patternfile_Count){PATTERNFILE=optarg;}
				else PATTERNFILE1=optarg;
				Patternfile_Count++;
				break;
			case 'o':
				JUNCTIONFILE=optarg;
				break;
			case 'w':
				WIGGLEFILE=optarg;
				break;
			case 'O':
				WRITE_SPLITREAD=TRUE;
				if (optarg) MAPFILE=optarg;
				break;
			case 'g':
				Name=optarg;Last_Dash=0;Genome_Name=optarg;Genome_String=optarg;
				for(;Name[0]!=0;Name++)
				{
					if (Name[0]=='/') 
					{
						Last_Dash++;Genome_Name=Name;
					}
				}

				REVBWTINDEX = (char*)Source;
				if(Last_Dash) Last_Dash=Genome_Name-optarg+1; else Genome_Name--;
				strncpy(REVBWTINDEX,optarg,Last_Dash);
				REVBWTINDEX[Last_Dash+0]='r';REVBWTINDEX[Last_Dash+1]='e';REVBWTINDEX[Last_Dash+2]='v';
				strcpy(REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(REVBWTINDEX+Last_Dash+3,".bwt"); 

				BWTFILE=REVBWTINDEX+500;
				strncpy(BWTFILE,optarg,Last_Dash);
				strcpy(BWTFILE+Last_Dash,Genome_Name+1);
				strcat(BWTFILE+Last_Dash,".bwt"); 


				REVOCCFILE = BWTFILE+500;
				strncpy(REVOCCFILE,optarg,Last_Dash);
				REVOCCFILE[Last_Dash+0]='r';REVOCCFILE[Last_Dash+1]='e';REVOCCFILE[Last_Dash+2]='v';
				strcpy(REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(REVOCCFILE+Last_Dash+3,".fmv"); 


				OCCFILE=REVOCCFILE+500;			
				strncpy(OCCFILE,optarg,Last_Dash);
				strcpy(OCCFILE+Last_Dash,Genome_Name+1);
				strcat(OCCFILE+Last_Dash,".fmv"); 

				SAFILE=OCCFILE+500;			
				strncpy(SAFILE,optarg,Last_Dash);
				strcpy(SAFILE+Last_Dash,Genome_Name+1);
				strcat(SAFILE+Last_Dash,".sa");

				REVSAFILE = SAFILE+500;
				strncpy(REVSAFILE,optarg,Last_Dash);
				REVSAFILE[Last_Dash+0]='r';REVSAFILE[Last_Dash+1]='e';REVSAFILE[Last_Dash+2]='v';
				strcpy(REVSAFILE+Last_Dash+3,Genome_Name+1);
				strcat(REVSAFILE+Last_Dash+3,".sa"); 

				BINFILE=REVSAFILE+500;			
				strncpy(BINFILE,optarg,Last_Dash);
				strcpy(BINFILE+Last_Dash,Genome_Name+1);
				strcat(BINFILE+Last_Dash,".pac");

				LOCATIONFILE=BINFILE+500;			
				strncpy(LOCATIONFILE,optarg,Last_Dash);
				strcpy(LOCATIONFILE+Last_Dash,Genome_Name+1);
				strcat(LOCATIONFILE+Last_Dash,".ann.location");

                                BLKFILE = LOCATIONFILE+500;
                                strncpy(BLKFILE,optarg,Last_Dash);
                                strcpy(BLKFILE+Last_Dash,Genome_Name+1);
                                strcat(BLKFILE+Last_Dash,".blk.");
				sprintf(BLKFILE+strlen(BLKFILE),"%d",RQFACTOR);

                                INDFILE = BLKFILE+500;
                                strncpy(INDFILE,optarg,Last_Dash);
                                strcpy(INDFILE+Last_Dash,Genome_Name+1);
                                strcat(INDFILE+Last_Dash,".ind.");
				sprintf(INDFILE+strlen(INDFILE),"%d",RQFACTOR);

                                RANGEFILE = INDFILE+500;
                                strncpy(RANGEFILE,optarg,Last_Dash);
                                strcpy(RANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(RANGEFILE+Last_Dash,".range");

                                SORTEDRANGEFILE = RANGEFILE+500;
                                strncpy(SORTEDRANGEFILE,optarg,Last_Dash);
                                strcpy(SORTEDRANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(SORTEDRANGEFILE+Last_Dash,".sort");
				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
	Patternfile_Count--;
	if (MARKEX) {Mark_Exons();exit(0);}
}

//}-----------------------------  Parse Command Line  -------------------------------------------------

//{-----------------------------  Misc  -------------------------------------------------
void Verbose()
{
	printf("\n-= RawBin Beta  =-\n");
	printf("Using the genome files\n %s\t %s\n %s\t %s\n", BWTFILE,OCCFILE,REVBWTINDEX,REVOCCFILE); 
	if (COMPRESS) printf ("Using a compressed index...\n");
	//printf("Query File : %s \t\t Output file: %s\n",PATTERNFILE,HITSFILE);
	printf("Query File : %s\n",PATTERNFILE);
	if (USEREFGENE) printf("Use refGene \n"); 

}

static const char LogTable256[] = 
{

  1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  
};

unsigned log2(unsigned v) // 32-bit word to find the log of
{
#ifdef BUILTIN_LOG
	return (32 -__builtin_clz(v));
#else
	unsigned r;     // r will be lg(v)
	register unsigned int t, tt; // temporaries

	if (tt = v >> 16)
	{
		r= (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	}
	else 
	{
		r= (t = v >> 8) ? 8 + LogTable256[t] : LogTable256[v];
	}
	return r;
#endif
}

//}-----------------------------  Misc  -------------------------------------------------

//{--------------------------------  Pair reads -------------------------------------------------------
void Pair_Reads(SARANGE *Head_Hits,SARANGE *Tail_Hits)
{
	SARANGE Head,Tail;
	int Pairs_Index,TGap,HGap;
	struct Valid_Hed
	{
		unsigned Location;
		bool Valid;
	}
	Multi_Head[SAGAP_CUTOFF];

	//for(int Orientation=0;Orientation<2;Orientation++)
	//{
	//if(Orientation) {Head_Hits=Head_Hits_Pos;Tail_Hits=Tail_Hits_Pos;} else {Head_Hits=Head_Hits_Neg;Tail_Hits=Tail_Hits_Neg;}
	for(int i=0;Head_Hits[i].Start;i++)//Iterate Head Hits
	{
		for(int j=0;Tail_Hits[j].Start;j++)//With Tail Hits
		{
			Head=Head_Hits[i];
			Tail=Tail_Hits[j];
			TGap=Tail.End-Tail.Start;
			HGap=Head.End-Head.Start;
			Pairs_Index=0;

			//if (TGap<= SAGAP_CUTOFF_INDEX || HGap <= SAGAP_CUTOFF)//Small sa ranges
			if (TGap<= SAGAP_CUTOFF_T || HGap <= SAGAP_CUTOFF_H)//Small sa ranges
			{
				if(TGap && HGap)//Head and tail multi hits...
				{
					unsigned Start;
					if(HGap<TGap)// || HGap > SAGAP_CUTOFF_INDEX)///debug=1
					{
						Start=Head.Start;
						for (int i=0;i<=HGap;i++)
						{
							Head.Start=Start;
							unsigned Loc=Head.End=Head.Start=Head.Conversion_Factor-BWTSaValue(revfmi,Head.Start);
							//Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							Multi_Head[i].Location=Loc;Multi_Head[i].Valid=true;
							for (int j=0;j<i;j++)
							{
								if (dist(Multi_Head[j].Location,Head.Start) <400000)
								{
									Multi_Head[j].Valid=false;Multi_Head[i].Valid=false;
								}
							}
							Start++;
						}
						for (int i=0;i<=HGap;i++)
						{
							if (Multi_Head[i].Valid)
							{
								Head.Start=Head.End=Multi_Head[i].Location;
								Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							}
						}
					}
					else
					{
						Start=Tail.Start;
						for (int i=0;i<=TGap;i++)
						{
							Tail.Start=Start;
							unsigned Loc=Tail.End=Tail.Start=Tail.Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
							//Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							Multi_Head[i].Location=Loc;Multi_Head[i].Valid=true;
							for (int j=0;j<i;j++)
							{
								if (dist(Multi_Head[j].Location,Tail.Start) <400000)
								{
									Multi_Head[j].Valid=false;Multi_Head[i].Valid=false;
								}
							}
							Start++;
						}
						for (int i=0;i<=TGap;i++)
						{
							if (Multi_Head[i].Valid)
							{
								Tail.Start=Tail.End=Multi_Head[i].Location;
								Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							}
						}
					}
				}
				else//Head or Tail unique..
				{
					Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
				}
			}
			else//Two big SA gaps...
			{
Sum++;T=1;
				//if (HGap > 500) continue;
//Sum++;T=1;
				//Get_Head_Tail(Head, Tail,EXONGAP,Pairs,Pairs_Index);
			}

			/*for(int i=0;Pairs[i].Head;i++)
			  {
			  Location_To_Genome(Pairs[i].Head);
			  char* ChrH=Genome_Offsets[Genome_Position].Genome;
			  Location_To_Genome(Pairs[i].Tail);
			  char* ChrT=Genome_Offsets[Genome_Position].Genome;
			  printf( "%u\t%s:%u\t%s\n",Pairs[i].Head,ChrH,Pairs[i].Tail,ChrT);
			  Total_Paired++;
			  }*/
		}
	}
	//}

}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Range_Index
 *  Description:  Loads the range index sturcture..
 * =====================================================================================
 */
void Load_Range_Index(char* INDFILE, char* BLKFILE,RANGEINDEX & Range_Index)
{
	unsigned Index_Size,Block_Size;

	Range_Index.Index=File_Open(INDFILE,"rb");
	Range_Index.Blocks=File_Open(BLKFILE,"rb");
	Range_Index.SA_Index=(SA*)malloc((Index_Size=Get_File_Size(Range_Index.Index)));//contains the index to blocks...
	Range_Index.SA_Blocks=(int*)malloc((Block_Size=Get_File_Size(Range_Index.Blocks)));
	if (!Range_Index.SA_Index || !Range_Index.SA_Blocks)
	{
		printf ("Load_Range_Index(): Memory allocation failed!..\n");
		exit(0);
	}
	fread(&Range_Index.COMPRESS,1,1,Range_Index.Blocks);
	fread(&Range_Index.Hits,1,sizeof(Range_Index.Hits),Range_Index.Blocks);
	fread(Range_Index.SA_Index,Index_Size,1,Range_Index.Index);
	fread(Range_Index.SA_Blocks,Block_Size,1,Range_Index.Blocks);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Small_Gap
 *  Description:  Do the searching when at least one of the pairs have small SA range... 
 * 		  Store the result in Pairs, starting from Pairs[Pairs_Index]
 * 		  Terminates hit list when encountering Pairs[x].Head=0... 
 * =====================================================================================
 */
void Search_Small_Gap(SARANGE & Head, SARANGE & Tail, unsigned d,PAIR* Pairs,int & Pairs_Index)
{
	unsigned H1,H2,T1,T2,L,H,M;
	TAG_INFO Head_Info,Tail_Info;


	if(Head.Start==Head.End)//Head is unique...
	{

		H1=Head.Start;//Conversion_Factor-BWTSaValue(revfmi,Head.Start);
		if (Tail.Start==Tail.End)// both hits are unique...
		{
			T2=Tail.Start;//Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
			if (T2>H1 && T2 < H1+d)//modify for multiple hits...
			{
				Pairs[Pairs_Index].Head=H1;
				Pairs[Pairs_Index].HLevel=Head.Level;
				Pairs[Pairs_Index].Tail=T2;
				Pairs[Pairs_Index].TLevel=Tail.Level;
				Pairs_Index++;
				if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
				if(H1>T2) 
				{
					printf("Search_Small_Gap(6):Enum error...\n");
					exit(0);
				}
#endif
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
		else //tail has multiples...
		{
			if (Tail.End-Tail.Start > SAGAP_CUTOFF)
			{
				Load_Info(Tail_Info,Tail);
				T1=Tail_Info.First;T2=Tail_Info.Last;


				if(H1<T1 )//Possible case for T1>H1 
				{
					if(H1+d>T1)
					{
						M=0;
						while (H1<T1 && H1+d>T1)//enumerate hits...
						{
							Pairs[Pairs_Index].Head=H1;
							Pairs[Pairs_Index].HLevel=Head.Level;
							Pairs[Pairs_Index].Tail=T1;
							Pairs[Pairs_Index].TLevel=Tail.Level;
							Pairs_Index++;
							if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(5):Enum error...\n");
							exit(0);
						}
#endif
							if (Pairs_Index>MAXCOUNT) 
							{
								Pairs[Pairs_Index].Head=0;
								return;
							}
							M++;
							T1=Get_Location(Tail_Info,M);
							if(M>=Tail_Info.Gap) break;
						}
					}
				}
				else if(T2>H1) //H1 inside tail gaps..
				{

					L=0;
					H=Tail_Info.Gap;
					while (L < H)
					{
						M=(L+H)/2;
						if (Get_Location(Tail_Info,M) > H1)
						{
							H=M;
						}
						else
						{
							L=M+1;
						}
					}
					if (L==H) M=H;//find tail position closest to unique head...

					T1=Get_Location(Tail_Info,M);
					while (H1<T1 && H1+d>T1)//enumerate hits...
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].HLevel=Head.Level;
						Pairs[Pairs_Index].Tail=T1;
						Pairs[Pairs_Index].TLevel=Tail.Level;
						Pairs_Index++;
						if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(4):Enum error...\n");
							exit(0);
						}
#endif
						if (Pairs_Index>MAXCOUNT) 
						{
							Pairs[Pairs_Index].Head=0;
							return;
						}
						M++;
						T1=Get_Location(Tail_Info,M);
						if(M>=Tail_Info.Gap) break;
					}
				}

				Pairs[Pairs_Index].Head=0;
				return;
			}
			else//Unique head and tail with gap below cutoff...
			{
				for (unsigned i=Tail.Start;i<=Tail.End;i++)
				{
					unsigned Hit=(Tail.Conversion_Factor-BWTSaValue(revfmi,i));
					if(Hit > H1 && Hit < H1+d)
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].HLevel=Head.Level;
						Pairs[Pairs_Index].Tail=Hit;
						Pairs[Pairs_Index].TLevel=Tail.Level;
						Pairs_Index++;
						if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>Hit) 
						{
							printf("Search_Small_Gap(3):Enum error...\n");
							exit(0);
						}
#endif

						if(Pairs_Index >MAXCOUNT) break;
					}
				}
				Pairs[Pairs_Index].Head=0;
				return;
			}
		}
	}
	else //if(Tail.End==Tail.Start)//Tail is unique...
	{
		T1=Tail.Start;//Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
		if(Head.End-Head.Start>SAGAP_CUTOFF)//Unique tail, but with multiple possible heads...
		{
			Load_Info(Head_Info,Head);
			H1=Head_Info.First;H2=Head_Info.Last;
			if(T1 > H1) //Head should not be after T1...
			{
				if(H2+d >T1)//Tail not too far away from heads...
				{
					//T1-d is between H1 and H2, search for the closest hit...
					L=0;
					H=Head_Info.Gap-1;
					T1-=d;
					while (L < H)
					{
						M=(L+H)/2;
						if (Get_Location(Head_Info,M) < T1)
						{
							L=M+1;
						}
						else
						{
							H=M;
						}
					}
					if (L==H) M=L;
					H1=Get_Location(Head_Info,M);

					T1+=d;
					while (H1<T1 && H1+d>T1)//enumerate hits...
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].HLevel=Head.Level;
						Pairs[Pairs_Index].Tail=T1;
						Pairs[Pairs_Index].TLevel=Tail.Level;
						Pairs_Index++;
						if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(2):Enum error...\n");
							exit(0);
						}
#endif
						if (Pairs_Index>MAXCOUNT) 
						{
							break;	
						}
						M++;
						H1=Get_Location(Head_Info,M);
						if(M>=Head_Info.Gap) break;
					}

				}
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
		else//Unique tail and head with gap below cutoff...
		{
			for (unsigned i=Head.Start;i<=Head.End;i++)//try all hits...
			{
				unsigned Hit=(Head.Conversion_Factor-BWTSaValue(revfmi,i));
				if(Hit < T1 && T1 < Hit+d)
				{
					Pairs[Pairs_Index].Head=Hit;
					Pairs[Pairs_Index].HLevel=Head.Level;
					Pairs[Pairs_Index].Tail=T1;
					Pairs[Pairs_Index].TLevel=Tail.Level;
					Pairs_Index++;
					if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
				if(Hit>T1) 
				{
					printf("Search_Small_Gap(1):Enum error...\n");
					exit(0);
				}
#endif
					if(Pairs_Index >MAXCOUNT) break;
				}
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Head_Tail...
 *  Description:  find the possible pairs for head and tail that are a distance d apart.. 
 *  		  Pairs are stored in the structure Pair..
 *  		  Pairs_Index is the offset of the sentinel indicating last hit...
 * =====================================================================================
 */

void Get_Head_Tail(SARANGE & Head, SARANGE & Tail,unsigned d,PAIR* Pairs,int & Pairs_Index)
{
	unsigned H1,H2,T1,T2;
	unsigned H,M,L;
	unsigned Location=0,Hd,Tl;

	TAG_INFO Head_Info,Tail_Info;


	Load_Info(Head_Info,Head);
	Load_Info(Tail_Info,Tail);
	H1=Head_Info.First;H2=Head_Info.Last + d;
	T1=Tail_Info.First;T2=Tail_Info.Last;


	if (T1>H2 || T2<H1) //disjoint head and pair...
	{
		Pairs[Pairs_Index].Head=0;//cannot be paired...
		return;
	}

	if (T1>H1) M=0;//Get in M the index of the closest value of Tail larger than H1
	else
	{
		L=0;
		H=Tail_Info.Gap-1;
		while (L < H)
		{
			M=(L+H)/2;
			if (Get_Location(Tail_Info,M) < H1)
			{
				L=M+1;
			}
			else
			{
				H=M;
			}
		}
		if (L==H) M=H;
	}


	Hd=H1;H=1;
	Tl=Get_Location(Tail_Info,M);M++;

#ifdef DEBUG
	if(Tl<Hd) {printf("Get_Head_Tail(): bin search error..\n");exit(0);}
#endif

	int MHead=1;//Array pos. of head...
	int TempM;
	unsigned TempTl;

	while (Hd<T2)//Tail must always be after Head...
	{
		if (Hd<Tl)
		{
			TempM=M;TempTl=Tl;
			while (Hd+d>Tl)//enum. Hits...
			{
				Pairs[Pairs_Index].Head=Hd;
				Pairs[Pairs_Index].HLevel=Head.Level;
				Pairs[Pairs_Index].Tail=Tl;
				Pairs[Pairs_Index].TLevel=Tail.Level;
				Pairs_Index++;
				if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
				if(Hd>Tl) 
				{
					printf("Get_Head_Tail():Enum error...\n");
					exit(0);
				}
#endif
				if(M >= Tail_Info.Gap) break;else Tl=Get_Location(Tail_Info,M);
				M++;
			}
			M=TempM;Tl=TempTl;

			if (MHead>=Head_Info.Gap) break; else Hd=Get_Location(Head_Info,MHead);
			MHead++;
		}
		else
		{
			if(M >= Tail_Info.Gap) break;else Tl=Get_Location(Tail_Info,M);
			M++;
		}
	}
	Pairs[Pairs_Index].Head=0;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Info
 *  Description:  Loads information necessary for reading coded info.... 
 * =====================================================================================
 */

void Load_Info( TAG_INFO & Tag, SARANGE & Head)
{
	Tag.Gap=Head.End-Head.Start+1;
	Tag.SA_Start=Head.Start;
	Tag.Block_Start=Get_Block_Start(Head.Start, Tag.Index);//get info block..
	Tag.First=SA_Index[Tag.Index].Start_Location;
	Tag.Last=SA_Index[Tag.Index].End_Location;
	Tag.Field_Length=log2(Tag.Gap);
}


inline unsigned Get_Location(TAG_INFO & Tag, unsigned Offset)
{
	unsigned SAPos;
	if(0==Offset) return Tag.First;
	if (Offset==Tag.Gap-1) return Tag.Last;
	Offset--;
	if(COMPRESS)
	{
		SAPos=Tag.SA_Start + bfx((unsigned char*)SA_Blocks,Tag.Block_Start+(Offset*Tag.Field_Length),Tag.Field_Length);
		return CONVERSION_FACTOR-BWTSaValue(revfmi,SAPos);
	}
	else
	{
		return (unsigned)SA_Blocks[Tag.Block_Start+Offset];
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Block_Start
 *  Description:  Finds the Block Value of the SAValue by bin search...
 *  		  Assumes SAValue exists...
 * =====================================================================================
 */
unsigned Get_Block_Start(unsigned SAValue,unsigned & M)
{
	unsigned L=0;
	unsigned H=Hits;

	while (L < H)
	{
		M=(L+H)/2;
		if (SA_Index[M].Start < SAValue)
		{
			L=M+1;
		}
		else
		{
			H=M;
		}
	}
	if (L==H) M=H; 

#ifdef DEBUG
	if (SA_Index[M].Start!=SAValue) 
	{
		printf("Get_Block_Start(): Bad SAValue !\n");
		exit(0);
	}
#endif

	return SA_Index[M].End;
}

//}--------------------------------  Pair reads -------------------------------------------------------

//{

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Location
 *  Description:  Loads genome locations file and returns the number of chromosomes..
 * =====================================================================================
 */
int Load_Locations(Offset_Record* Genome_Offsets,unsigned Offsets[])
{
	FILE* Location_File=File_Open(LOCATIONFILE,"r");
	int Genome_Count=0;

	while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<80)
	{
		Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
		fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File);
		strcpy(Genome_Offsets[Genome_Count].GenomeM,Genome_Offsets[Genome_Count].Genome);
		for(int i=0;i<40;i++) 
		{
			if (Genome_Offsets[Genome_Count].Genome[i] == '\n' ||Genome_Offsets[Genome_Count].Genome[i] == '\r')
			{ 
				Genome_Offsets[Genome_Count].Genome[i]=0;
				Genome_Offsets[Genome_Count].GenomeM[i]='-';
				Genome_Offsets[Genome_Count].GenomeM[i+1]=0;
				break;
			} 
		}
		//if(MAPMODE) Genome_Offsets[Genome_Count].Out_File=File_Open(Genome_Offsets[Genome_Count].Genome,"w+b");
		//else Genome_Offsets[Genome_Count].Out_File=File_Open(Genome_Offsets[Genome_Count].Genome,"r");
		Genome_Count++;	
	}
	for ( int i=1;i<Genome_Count;i++)
	{
		Offsets[i]=Offsets[i-1]+Genome_Offsets[i].Offset;
		if(MAPMODE) 
		{
			string S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_File=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_FileM=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].Unmapped=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].UnmappedM=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedX=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedXM=File_Open(S.c_str(),"w+b");
			
		}
		else 
		{
			string S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_File=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_FileM=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].Unmapped=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].UnmappedM=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedX=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedXM=File_Open(S.c_str(),"r");
		}
		if (USEREFGENE) 
		{
			string S=Genome_String+Genome_Offsets[i-1].Genome+".ref";
			Genome_Offsets[i-1].Ref_File=File_Open(S.c_str(),"r");
		}
	}
	return Genome_Count-1;
}
//}

int Location_To_Genome(unsigned & Location)
{
	int Genome_Position=0;
	while ( Genome_Position< Genome_Count )
	{
		if (Location < Offsets[Genome_Position]) break;
		Genome_Position++;
	}
	Genome_Position--;
	Location=Location-Offsets[Genome_Position];
	return Genome_Position;
}

void Show_Progress(unsigned Percentage)
{
	if (Percentage >=97) return;
	printf("+%u%%\b\b\b",Percentage);
	fflush(stdout);
}

void Load_RefGene(FILE* Ref_Handle,Hash & Ref_Gene,JUNC_HASH & Donor_AT,JUNC_HASH & Donor_GT,JUNC_HASH & Donor_GC,JUNC_HASH & Donor_CT,JUNC_HASH & Acc_AC,JUNC_HASH & Acc_AG,JUNC_HASH & Acc_AT,JUNC_HASH & Acc_GC)
{
	char* String_Buffer=new char[MAX_REF_LINE+1];
	char** Fields=new char*[20];
	char** X_Start=new char*[MAX_X];
	char** X_End=new char*[MAX_X];
	char* & Chr=Fields[2];
	char* & Exon_Count=Fields[8];
	char* & Exon_Start=Fields[9];
	char* & Exon_End=Fields[10];
	char* & Strand=Fields[3];
	int Max_Ex=0;


	try
	{
		while(1)
		{
			if (fgets(String_Buffer,MAX_REF_LINE,Ref_Handle))
			{
				int Tab_Count=0;
				Split(String_Buffer,'\t',Fields);
				int Exon_Num=atoi(Exon_Count);
				Split(Exon_Start,',',X_Start,MAX_X);
				Split(Exon_End,',',X_End,MAX_X);
				if (Exon_Num > MAX_X) {printf("Long splicing ... Some exons skipped..\n");Exon_Num=MAX_X;}
				for (int i=0;i<Exon_Num;i++) 
				{
					//printf("%s\n",Strand);
					int Pairing=0;
					unsigned x=atoi(X_End[i])-1,y=atoi(X_Start[i]);
					if(*Strand=='+') //Pairing=De_Novo_Splice(JPair,JStat);
					{

						unsigned char AccA= Bit2Char(y-2);
						unsigned char AccB= Bit2Char(y-1);
						if (AccA == N_a)
						{
							if(AccB == N_g)//AG -high chance case [GT-AG] or [GC - AG]
							{ 
								unsigned char A,B,C,D;
								A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
								Acc_AG[y]= (D | (C<<2) | (B<<4) |(A<<6));
							}
							else if (AccB == N_c)//AC - [AT-AC]
							{
								unsigned char A,B,C,D;
								A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
								Acc_AC[y]= (D | (C<<2) | (B<<4) |(A<<6));
							}
						}
						unsigned char DonorA= Bit2Char(x+1);
						unsigned char DonorB= Bit2Char(x+2);
						if (DonorA == N_g)
						{
							if(DonorB == N_t)//GT [GT-AG]
							{
								unsigned char A,B,C,D;
								A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
								Donor_GT[x]= (D | (C<<2) | (B<<4) |(A<<6));
							}
							else if(DonorB == N_c)//GC [GC-AG]
							{
								unsigned char A,B,C,D;
								A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
								Donor_GC[x]= (D | (C<<2) | (B<<4) |(A<<6));
							}
						} 
						if (DonorA == N_a && DonorB == N_t)//AT [AT-AC]
						{
							unsigned char A,B,C,D;
							A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
							Donor_AT[x]= (D | (C<<2) | (B<<4) |(A<<6));
						}

					}
					else //Pairing=De_Novo_Splice_Rev(JPair,JStat); 
					{
						unsigned char AccA= Bit2Char(y-2);
						unsigned char AccB= Bit2Char(y-1);
						if (AccA == N_a)
						{
							if(AccB == N_c)//GC [CT-GC]
							{ 
								unsigned char A,B,C,D;
								A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
								Acc_AC[y]= (D | (C<<2) | (B<<4) |(A<<6));
							}
							else if (AccB == N_t)//GT - [AT-GT]
							{
								unsigned char A,B,C,D;
								A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
								Acc_AT[y]= (D | (C<<2) | (B<<4) |(A<<6));
							}
						}
						if (AccA == N_g && AccB == N_c)//AC [CT-AC]
						{
							unsigned char A,B,C,D;
							A=Bit2Char(y);B=Bit2Char(y+1);C=Bit2Char(y+2);D=Bit2Char(y+3);
							Acc_GC[y]= (D | (C<<2) | (B<<4) |(A<<6));
						}
						//------------------------------ Seek canonical Donor sites ------------------------------------------------------
						unsigned char DonorA= Bit2Char(x+1);
						unsigned char DonorB= Bit2Char(x+2);
						if (DonorB == N_t)
						{
							if(DonorA == N_c)//CT [CT-AC][CT -GC]
							{
								unsigned char A,B,C,D;
								A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
								Donor_CT[x]= (D | (C<<2) | (B<<4) |(A<<6));
								Pairing= S_ct;
							}
							else if(DonorA == N_g)//GT [AT-GT]
							{
								unsigned char A,B,C,D;
								A=Bit2Char(x-3);B=Bit2Char(x-2);C=Bit2Char(x-1);D=Bit2Char(x);
								Donor_GT[x]= (D | (C<<2) | (B<<4) |(A<<6));
								Pairing= S_gt;
							}
						} 

					}
				}

			}
			else if(!feof(Ref_Handle)) throw("Read Error..\n"); else break;
		}
	}
	catch(char *Error)
	{
		printf("%s",Error);
		exit(0);
	}
	delete [] String_Buffer;delete [] Fields;delete [] X_End;delete []X_Start;
	printf("%d\n",Max_Ex);

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Split
 *  Description:  split string into a maximum of Max fields seperated by Sep, and strore starts in Fields..
 *  		  if Max=0, all fields are split..
 *  		  returns the number of fields..
 * =====================================================================================
 */
int Split(char* String,char Sep, char* Fields[],int Max)
{
	int i=0,j=0;
	char* Last=String;
	while(String[i] && String[i]!='\n')
	{
		if(String[i]==Sep)
		{
			Fields[j++]=Last;
			String[i]=0;Last=String+i+1;
			if (Max ==j) break;
		}
		i++;
	}
	return j;
}


bool Scan(unsigned & Total_Hits,unsigned & Tags_Processed,int Off,unsigned Read_ID)
{
	LH=STRINGLENGTH/2;//calculate tag portions...
	if ((STRINGLENGTH % 2)) {LH++;}	
	RH=STRINGLENGTH-LH;
	RHQL=RH/2;RHQR=RH-RHQL;

	Guess_Orientation(Off);

//if (HITS) return true;
	//if (HITS==1 && SKIP_INIT_HIT) {Total_Hits++;return true;}//make this skip on unit hit?
	//{----------------------- Scan half and extend --------------------------------------------------
	char Hit_Found=FALSE;
	char Strand=0;
	char L=0,R=1;

	while (!Hit_Found && Strand<1)//at most two exons possible?
	{
		if (Orientation[Strand].Strand=='-') {L=1;R=0;}
		{
			if (Try[LEFT])//left half has an exact match..
			{
				if(Scan_Half_ExtendR(Orientation[L].String,Orientation[L],Read_ID))//try exon in right half..
				{Hit_Found=TRUE;Tags_Processed++;}//continue;}
			}
			if (Try[RIGHT])
			{
				if(Scan_Half_ExtendL(Orientation[L].String,Orientation[L],Read_ID)) 
				{Hit_Found=TRUE;Tags_Processed++;}//continue;}
			}
			if (Try[LEFTC])//left half has an exact match..
			{
				if(Scan_Half_ExtendR(Orientation[R].String,Orientation[R],Read_ID))//try exon in right half..
				{Hit_Found=TRUE;Tags_Processed++;}//continue;}
			}
			if (Try[RIGHTC])
			{
				if(Scan_Half_ExtendL(Orientation[R].String,Orientation[R],Read_ID)) 
				{Hit_Found=TRUE;Tags_Processed++;}//continue;}
			}
			//mishits ... printf("*\t%s",Head.Description+1);
			Strand++;
			continue;// no two exons...
		}
	}
	return Hit_Found;
}

void Write_Sam(char Sign,char* Chr,unsigned Loc1,unsigned Loc2,int G1,int G2,char Code,unsigned Read_ID)
{
	unsigned char Flag=0;
	char Q[MAXSTRINGLEN];
	char C[MAXSTRINGLEN];
	char* QString=Head.Quality;
	char* Seq= Head.Tag_Copy;
	char* CIGAR;
	//assert(Loc1 <=Loc2 || (Loc2==0 && G2==0));
	READ HeadT=Head;

	if (Sign == '-') 
	{
		Flag =16;
		Seq=C;
		for (int i=0;i<STRINGLENGTH;i++) C[i]=Code_To_CharCAPS[Head.Complement[i]];C[STRINGLENGTH]=0;
		if(FILETYPE == FQ) {for (int i=0;i<STRINGLENGTH;i++) Q[i]=Head.Quality[STRINGLENGTH-i];QString=Q;}

	}
	char* S;for (S=Head.Description;*S;S++);//search for string end..
	char T1=*(S-1);*(S-1)=0;//remove CR
	char T2=Seq[STRINGLENGTH];Seq[STRINGLENGTH]=0;
	if(G2)//split map
	{
		fprintf (SAM,"%s\t%d\t%s\t%u\t255\t%dM%dN%dM\t*\t0\t0\t%s\t%s\tXT:i:%c\tXL:i:%d\tXI:i:%u\n",Head.Description+1,Flag,Chr,Loc1,G1,Loc2-(Loc1+G1),G2,Seq,QString,Code,STRINGLENGTH,Read_ID);//wrong CIGAR
	}
	else
	{
		fprintf (SAM,"%s\t%d\t%s\t%u\t255\t%dM\t*\t0\t0\t%s\t%s\tXT:i:%c\tXL:i:%d\tXI:i:%u\n",Head.Description+1,Flag,Chr,Loc1,G1,Seq,QString,Code,STRINGLENGTH,Read_ID);
	}
	*(S-1)=T1;//remove CR
	Seq[STRINGLENGTH]=T2;
	//assert(HeadT==Head);
	Head=HeadT;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Juncs_In_Read(int Head_Start,int Head_Length,int Tail_Start, int Tail_Length)
 *  Description:  See if there are junction motifs in the read, and see if they are more than the allowed number..
 * =====================================================================================
 */
int Juncs_In_Read(OPX* Junc_List,int & Junc_List_Ptr,int Head_Start,int Head_Length,int Tail_Start, int Tail_Length)
{
	assert (Junc_List && Junc_List_Ptr < MAX_JUNC_LIST && Head_Start >0 && Head_Length >0 && Tail_Start >0 && Tail_Length >0);
	int Possible_Junc_Count=0;
	bool Motif_Found=false;
	if(Head_Length && Tail_Length && Head_Length + Tail_Length >= STRINGLENGTH)
	{
		int Gap=Head_Length-(STRINGLENGTH-Tail_Length);
		int Motif=0;
		assert(Gap>=0);//if (Gap >4) return 0;
		for (int i=0;i<=Gap;i++) 
		{
			if(Motif=Scan_Plus_Motif(Head_Start+STRINGLENGTH-Tail_Length+i-1,Tail_Start+i))
			{
				Junc_List[Junc_List_Ptr].x=Head_Start+STRINGLENGTH-Tail_Length+i-1;
				Junc_List[Junc_List_Ptr].y=Tail_Start+i;
				Junc_List[Junc_List_Ptr++].Motif=Motif;
				Motif_Found=true;
				if(++Possible_Junc_Count > MAX_JUNC_IN_READ) break;
			}
			else if(Motif=Scan_Minus_Motif(Head_Start+STRINGLENGTH-Tail_Length+i-1,Tail_Start+i))
			{
				Junc_List[Junc_List_Ptr].x=Head_Start+STRINGLENGTH-Tail_Length+i-1;
				Junc_List[Junc_List_Ptr].y=Tail_Start+i;
				Junc_List[Junc_List_Ptr++].Motif=Motif;
				Motif_Found=true;
				if(++Possible_Junc_Count > MAX_JUNC_IN_READ) break;
			}
		}
	}
	if (Junc_List_Ptr >= MAX_JUNC_LIST) Junc_List_Ptr=MAX_JUNC_LIST; //This does not increment the Junc_List_Ptr, which make this statement useless.
	return Possible_Junc_Count;
}


