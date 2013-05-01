#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
extern "C" 
{
	//#include "iniparser.h"
	//#include <time.h>
	//#include "MemManager.h"
	//#include "MiscUtilities.h"
	//include "TextConverter.h"
	#include "BWT.h"
}

using namespace std;

struct SARANGE
{
	/*unsigned Start;
	unsigned End;
	int Level;//at what relative node are we?
	char Mismatches;//number of mismatches at this level
	unsigned char Skip;
	unsigned Mismatch_Char;//2|2|...
	unsigned Conversion_Factor;*/
	

	unsigned Start;//is UINT_MAX if invalid hit
	unsigned End;
	int Level;//at what relative node are we?
	char Mismatches;//number of mismatches at this level
	unsigned char Skip;
};

const int MAXHITS=10;
const int SAINTERVAL=8;
const int BRANCHTHRESHOLD=80;//30 //Threshold at which to check the BWT instead of branching
const int MAX_EXON_LIMIT=18;
BWT *revfmi;
int RESIDUE=5;
int STRINGLENGTH;
unsigned SOURCELENGTH;
char Char_To_Code[256];
char Code_To_Complement[255];

//bool Search_Forwards_Exact(SARANGE & Tag,char* Current_Tag, int Start,int StringLength);//,unsigned Read_ID);//,SARANGE & Ret_SA);
void Search_Forwards_Exact(char *Current_Tag,SARANGE & Tag,int Start,int StringLength,BWT *fmi);
void Get_SARange( char New_Char,SARANGE & Range,BWT *fmi);
bool Call_Junc(string & L,string & DesS);
void Init_Char_To_Code(char *Char_To_Code);
void Bin_Encode(char *Current_Tag,char *Char_To_Code,int STRINGLENGTH);
void Load_Indexes(BWT* & fwfmi,unsigned & SOURCELENGTH,char *BWT,char *OCC,char *SA);
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName,const char* SAFile);
void Rev_Comp(char *Read,int StringLength);
void Dump_Fasta();
//>seq.1:chr18:82754660-82754735/1        CCATCACCTCCCTTCCTAGGCGTGACGGTGCCAGTGTGGGGCTCTTTGGAAGGTCTACAGAAGGCATGCTAAATCA 
main (int argc,char* argv[])
{
	
	char REVBWTFILE[]="/home/Ken/rawbin/index/batman/mm9/revmm9.bwt" ; 
	char REVOCCFILE[]= "/home/Ken/rawbin/index/batman/mm9/revmm9.fmv";
	//char REVSAFILE[]= "/home/Ken/rawbin/index/batman/mm9/revmm9.sa";

	if(argc<2)
	{
		printf("%s <batman_output> \n",argv[0]);
		exit(100);
	}
	try
	{
		//Load_Indexes(revfmi,SOURCELENGTH,REVBWTFILE,REVOCCFILE,NULL);//REVSAFILE);
		Init_Char_To_Code(Char_To_Code);
		for (int i=0;i<255;i++) Code_To_Complement[i]=0;for (int i=0;i<4;i++) Code_To_Complement[i]=3-i;

		string L,HitsE[MAXHITS],HitsB[MAXHITS],Des;
		bool In_Read=true;
		ifstream IN(argv[1]); if(!IN.is_open())throw ("Cannot open file ..."); 
		FILE *MISHIT;MISHIT=fopen("Multi_Hit.fa","w"); if (!MISHIT) throw ("Cannot create file...");
		getline(IN,L);
		if (L[0]!='@') throw ("Corrupt outut?..\n");

		while (IN.good())
		{
			int BHits=0,EHits=0,HitsEP=0,HitsBP=0;
			getline(IN,Des);
			In_Read=true;
			while(In_Read)
			{
				getline(IN,L);
				if(L[0]=='@')
				{
					In_Read=false;break;
				}
				else
				{
					if (L[0]=='*')
					{
						BHits++;
						HitsB[HitsBP++]=L;
					}
					else
					{
						EHits++;
						HitsE[HitsEP++]=L;
					}
					assert (HitsEP<MAXHITS && HitsBP<MAXHITS);
				}
				if(IN.eof()) {In_Read=false;}//File_Full=false;}
			}

			if(BHits)
			{
				if((BHits+EHits)==1)
				{
					if(Call_Junc(HitsB[0],Des))
					{
					}
				}
				else
				{
					char Description[1000],Read[500];
					sscanf(Des.c_str(),"%s%s",Description,Read);
					//MISHIT << Description << "\n" << Read <<"\n";
					fprintf(MISHIT,">%s\n%s\n",Description,Read);
				}
			}
		}
	}
	catch(const char* Err)
	{
		printf ("%s\n",Err);
	}
}

bool Call_Junc(string & L,string & DesS)
{
	assert (L.c_str() && DesS.c_str());
	char Des[1000],Junc[1000],Dummy[5000],Right[500],Left[500],Chr[20],Strand[3],Read[500],SignM;
	unsigned Pos,RightL,LeftL,RightR,LeftR;
	static bool First_Pass=true;

	sscanf(DesS.c_str(),"%s%s",Des,Read);
	if (First_Pass)
	{
		STRINGLENGTH=strlen(Read);
		First_Pass=false;
	}
	sscanf(L.c_str(),"%s%s%s%u",Dummy,Junc,Strand,&Pos);
	SignM=Strand[0];

	char *token = strtok (Junc, ",");
	sscanf (token, "%s", Left);
	token = strtok (NULL,",");
	if(token)
	{
		    sscanf (token, "%s", Right);
	}
	else {Right[0]=0;}

	if(!Right[0] || 'E'==Right[0]|| '<'==Right[0])
	{
		return 0;
	}
	else
	{
//*       chrX:7176283:7176462:-,chrX:7176559:7176626:-   -       105     0       76      692

		sscanf(Left,"%[^:]%*c%u%*c%u%s",Chr,&LeftL,&LeftR,Strand);
		sscanf(Right,"%*[^:]%*c%u%*c%u%s",&RightL,&RightR,Strand);
		char Sign=Strand[1];
		assert(LeftL>0 && LeftR >0 && RightL >0 && RightR>0 && (Sign=='+' || Sign=='-'));

		int Left_Gap=LeftR-(LeftL+Pos);
		assert(Left_Gap >0);
		if ((Left_Gap<RESIDUE) || (Left_Gap>STRINGLENGTH-RESIDUE))
		{
			return 0;
		}
		else
		{
			Bin_Encode(Read,Char_To_Code,STRINGLENGTH);
			if('-'==SignM) Rev_Comp(Read,STRINGLENGTH);

			SARANGE Range;
			bool Skip_Hit=false;
			/*if (Left_Gap >= STRINGLENGTH/2)
			{
				Range.Start=0;Range.End=SOURCELENGTH;Range.Level=1;Range.Mismatches=0;
				Search_Forwards_Exact(Read,Range,1,Left_Gap,revfmi);
				if(Range.Start != UINT_MAX)
				{
					if(Range.End-Range.Start >2) Skip_Hit=true;
				}	
				if(!Skip_Hit)
				{
					Rev_Comp(Read,Left_Gap);
					Range.Start=0;Range.End=SOURCELENGTH;Range.Level=1;Range.Mismatches=0;
					Search_Forwards_Exact(Read,Range,1,Left_Gap,revfmi);
					if(Range.Start != UINT_MAX)
					{
						if(Range.End-Range.Start >2) Skip_Hit=true;
					}	
				}
				else
				{
					Dump_Fasta();
				}
			}
			else
			{
				Range.Start=0;Range.End=SOURCELENGTH;Range.Level=1;Range.Mismatches=0;
				Search_Forwards_Exact(Read+Left_Gap,Range,1,STRINGLENGTH-Left_Gap,revfmi);
				if(Range.Start != UINT_MAX)
				{
					if(Range.End-Range.Start >2) Skip_Hit=true;
				}	
				if(!Skip_Hit)
				{
					Rev_Comp(Read+Left_Gap,STRINGLENGTH-Left_Gap);
					Range.Start=0;Range.End=SOURCELENGTH;Range.Level=1;Range.Mismatches=0;
					Search_Forwards_Exact(Read+Left_Gap,Range,1,STRINGLENGTH-Left_Gap,revfmi);
					if(Range.Start != UINT_MAX)
					{
						if(Range.End-Range.Start >2) Skip_Hit=true;
					}	
				}
				else
				{
					Dump_Fasta();
				}
			}*/
			
			if(!Skip_Hit)
			{
				int Left_Cord=LeftR-3;int Right_Cord=RightL+3;
				int Gap=Right_Cord-Left_Cord-3;
				assert(Gap >0);
				printf("%s\t%u\t%u\t%s\t1000\t%c\t%u\t%u\t255,0,\t2\t3,3\t0,%d\n",Chr,Left_Cord,Right_Cord,Des,Sign,Left_Cord,Right_Cord-1,Gap);
				return 1;
			}
			else return 0;
		}
	}
	return true;
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
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Forwards_Exact
 *  Description:  forward seach for exact occurence for string in Current_Tag in revfmi
 *  		  Start is the 1-indexed location of string position in Current_Tag
 *  		  String length is the length of substring to search...
 *  		  Start+Tag.Length-2 last position...
 *  		  return true if full StringLength scanned.
 * =====================================================================================
 */

bool Search_Forwards_Exact(SARANGE & Tag,char* Current_Tag, int Start,int StringLength)//,unsigned Read_ID)//,SARANGE & Ret_SA)
{
	unsigned Index,First,Last;
	unsigned Branch_Characters[4];
	SARANGE Temp;//temporary tag to save last tag details before failure...
	char Now;
	Tag.Level=1;
	Start -= 2;//accomodate 1 based offset.. Start is 1 based and so is tag Level
	//Ret_SA.Start=0;
	//Tag=Rev_Fmi[Current_Tag[Start+Tag.Level]];
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

					//Cache_SF[Tag.Level]=Tag;
					if(Tag.Level== StringLength)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						//if (!Ret_SA.Start) Ret_SA=Tag;
						return true; 	
					}
					else {Tag.Level++;continue;}
				} 
				else//mismatch...
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level--;
					//if (!Ret_SA.Start) Ret_SA=Tag;
					return false;	
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
					//if (!Ret_SA.Start && Tag.End-Tag.Start <= MAX_HITS_IN_SMALL_RES) Ret_SA=Tag;
				}
				else//mismatch..
				{
					//if (!Ret_SA.Start) Ret_SA=Tag;
					return false;
				}
			} 
			else
			{
				Temp=Tag;
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,revfmi);
				if (!Tag.Start) 
				{
					//Tag=Temp;
					//if (!Ret_SA.Start) Ret_SA=Tag;
					return false;
				}
				//if (!Ret_SA.Start && Tag.End-Tag.Start <= MAX_HITS_IN_SMALL_RES) Ret_SA=Tag;
			}

			//Cache_SF[Tag.Level]=Tag;
			if(Tag.Level== StringLength)
			{
				//if (!Ret_SA.Start) Ret_SA=Tag;
				return true;
			}
			else {Tag.Level++;continue;}

		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Forwards_Exact
 *  Description:  gets the SA* range of Current_Tag[Start+1...StringLength]; 
 *  Require:	  Tag=SARange*[1..Start] 
 * =====================================================================================
 */
void Search_Forwards_Exact(char *Current_Tag,SARANGE & Tag,int Start,int StringLength,BWT *fmi)
{
	Start-=2;
	assert(Tag.Start!=UINT_MAX);
	for(;;)	
	{
		Get_SARange(Current_Tag[Start+Tag.Level],Tag,fmi);

		if (Tag.Start!=UINT_MAX)
		{
			if(Tag.Level== StringLength)
			{
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else//Mismatch
		{
			return;
		}

	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initFMI
 *  Description:  Opens FM index fmiFile
 * =====================================================================================
 */
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName,const char* SAFile) 
{
	BWT *fmi;
	MMPool *mmPool;
        int PoolSize = 524288;
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	fmi = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SAFile, NULL, NULL, NULL);//Load FM index
	return fmi;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Indexes
 *  Description:  Load reverse and forward indexes...
 * =====================================================================================
 */
void Load_Indexes(BWT* & fwfmi,unsigned & SOURCELENGTH,char *BWT,char *OCC,char *SA)
{
	fwfmi=initFMI(BWT,OCC,SA);//Load FM indexes
	SOURCELENGTH = fwfmi->textLength;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Bin_Encode
 *  Description:  Convert Nucleotide string to binary format
 *  Require:	  Current_Tag-string consisting of Nuces.., Char_To_Code- Translator array,
 *  		  STRINGLENGTH-length of string ..
 * =====================================================================================
 */
void Bin_Encode(char *Current_Tag,char *Char_To_Code,int STRINGLENGTH)
{
	for (unsigned i=0;i<=STRINGLENGTH-1;i++)
	{
		Current_Tag[i]=Char_To_Code[Current_Tag[i]];
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Init_Char_To_Code
 *  Description:  Initialise nuce to bin translator
 *  Require:	  Char_To_Code - Converter array..
 * =====================================================================================
 */

void Init_Char_To_Code(char *Char_To_Code)
{
	for (int i=0;i<256;i++) Char_To_Code[i]=0;
	Char_To_Code['N']=0;Char_To_Code['n']=0;Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_SARange
 *  Description:  gets the SA range of strings having prefix [New_Char][Range]
 * =====================================================================================
 */

void Get_SARange( char New_Char,SARANGE & Range,BWT *fmi)
{
	assert(Range.Start!=UINT_MAX);
	Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.Start, New_Char) + 1;
	Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.End+1, New_Char);
	if (Range.End<Range.Start) 
	{
		Range.Start=UINT_MAX;
	}

}


void Rev_Comp(char *Read,int StringLength)
{
	char Complement[500];
	for (unsigned i=0;i<=StringLength-1;i++)
	{
		Complement[StringLength-1-i]=Code_To_Complement[Read[i]];//change later...
	}
	for (unsigned i=0;i<=StringLength-1;i++) Read[i]=Complement[i];
}

void Dump_Fasta()
{
	return;
}
