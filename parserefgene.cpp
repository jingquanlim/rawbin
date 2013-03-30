#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <string.h>
#include "ordpair.h"
using namespace std;
const int MAX_REF_LINE=10000;
const int MAX_X=500;//for mouse 312
struct exon
{
	unsigned Start,End;
};

struct Offset_Record
{
	char Genome[40];
	unsigned Offset;
	FILE* Out_File;
};

FILE* File_Open(const char* File_Name,const char* Mode);
int Split(char* String,char Sep, char* Fields[],int Max=0);
void Load_RefGene(char* File_Name,Hash & Ref_Gene);
int Filter_RefGene(char* Genome,char* File_Name);


int main( int argc,char* argv[])
{
	if (argc <3)
	{
		printf("Processes RefGene for Rawbin..\n");
		printf("%s <genome> <refgene>\n",argv[0]);
		exit(0);
	}
	Filter_RefGene(argv[1],argv[2]);
	//Load_RefGene(argv[1],Ref_Gene);
}

void Load_RefGene(char* File_Name,Hash & Ref_Gene)
{
	char* String_Buffer=new char[MAX_REF_LINE+1];
	char** Fields=new char*[20];
	char** X_Start=new char*[MAX_X];
	char** X_End=new char*[MAX_X];
	FILE* Ref_Handle=File_Open(File_Name,"r");
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
					printf("%s\n",Strand);
					//OP.Start=X_Start[i];OP.End=X_End[i];
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


int Filter_RefGene(char* Genome,char* File_Name)
{
	int Genome_Count=0;
	Offset_Record Genome_Offsets[80];
	map <string,FILE*> GList;
	map <string,FILE*> ::iterator GList_IT;
	char* String_Buffer=new char[MAX_REF_LINE+1];
	char* String_Org=new char[MAX_REF_LINE+1];
	char** Fields=new char*[20];
	char* & Chr=Fields[2];
	FILE* Ref_Handle=File_Open(File_Name,"r");
	string LOCATIONFILE=Genome;
	LOCATIONFILE+=".ann.location";
	FILE* Location_File=File_Open(LOCATIONFILE.c_str(),"r");

	while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<80)
	{
		Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
		fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File);
		strcpy(Genome_Offsets[Genome_Count].Genome,Genome_Offsets[Genome_Count].Genome);
		for(int i=0;i<40;i++) 
		{
			if (Genome_Offsets[Genome_Count].Genome[i] == '\n' ||Genome_Offsets[Genome_Count].Genome[i] == '\r')
			{ 
				Genome_Offsets[Genome_Count].Genome[i]=0;
				break;
			} 
		}
		Genome_Count++;	
	}
	for ( int i=1;i<Genome_Count;i++)
	{
		string S=Genome_Offsets[i-1].Genome;
		S=Genome_Offsets[i-1].Genome;S=Genome+S+".ref";
		Genome_Offsets[i-1].Out_File=File_Open(S.c_str(),"w");
		GList[Genome_Offsets[i-1].Genome]=Genome_Offsets[i-1].Out_File;
	}

	try
	{
		while(1)
		{
			if (fgets(String_Buffer,MAX_REF_LINE,Ref_Handle))
			{
				strcpy(String_Org,String_Buffer);
				int Tab_Count=0;
				Split(String_Buffer,'\t',Fields);
				GList_IT=GList.find(Chr);
				if(GList_IT != GList.end())
				{
					fprintf(GList_IT->second,"%s",String_Org);
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
	delete [] String_Buffer;delete [] Fields;
	return Genome_Count-1;
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

char* Get_Field(char* String_Buffer,int & Tab_Count,int Skip)
{
return NULL;
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
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}
