#include "zlib.h"
#include "stdio.h"
#include "stdlib.h"
#include "common.h"
#include "const.h"
#include "file.h"
#include "assert.h"

extern char Char_To_CodeC[];
extern char Char_To_Code[];
extern int READLEN;
extern int ORG_STRINGLENGTH;
extern bool SAM_READER;
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

FILE* File_Exist_Open(const char* File_Name)
{
	FILE* Handle;
	Handle=fopen64(File_Name,"r");
	return Handle;
}

unsigned Get_File_Size(FILE* File)
{
	fseek (File , 0 , SEEK_END);
	unsigned Size = ftell (File);
	rewind (File);
	return Size;
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
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Detect_Input
 *  Description:  Detects file type and lengths of input file(s)
 * =====================================================================================
 */
void Detect_Input(FILETYPE & P, FILE* & Input_File,FILE* & Mate_File)
{
	char Description[MAXDES+1];
	char Current_Tag[MAXTAG+1];
	char Quality[MAXTAG+1];

	if (fgets(Description,MAXDES,Input_File)!=0)//Measure tag length
	{
		fgets(Current_Tag,MAXTAG,Input_File);
		for(P.TAG_COPY_LEN=0;Current_Tag[P.TAG_COPY_LEN]!='\n' && Current_Tag[P.TAG_COPY_LEN]!='\r' && Current_Tag[P.TAG_COPY_LEN]!=0;P.TAG_COPY_LEN++);//TAG_COPY_LEN++;
		if(Patternfile_Count)
		{
			fgets(Description,MAXDES,Mate_File);
			Current_Tag[P.TAG_COPY_LEN++]='\t';fgets(Current_Tag+P.TAG_COPY_LEN,MAXTAG,Mate_File);
			for(P.TAG_COPY_LEN=0;Current_Tag[P.TAG_COPY_LEN]!='\n' && Current_Tag[P.TAG_COPY_LEN]!='\r' && Current_Tag[P.TAG_COPY_LEN]!=0;P.TAG_COPY_LEN++);//TAG_COPY_LEN++;
		}
		for(P.STRINGLENGTH=0;Current_Tag[P.STRINGLENGTH]!='\n' && Current_Tag[P.STRINGLENGTH]!='\r' && Current_Tag[P.STRINGLENGTH]!=0 && Current_Tag[P.STRINGLENGTH]!=PAIR_END_SEPERATOR;P.STRINGLENGTH++);
		if(Current_Tag[P.STRINGLENGTH]==PAIR_END_SEPERATOR) 
		{
			P.NORMAL_TAGS=FALSE;//we have pair ended tags..
			if(Patternfile_Count) {P.PAIRING_TYPE=TWOFILE;}else {P.PAIRING_TYPE=TAB;}
			P.PAIR_LENGTH_LEFT=P.STRINGLENGTH;

			for(P.PAIR_LENGTH_RIGHT=0;Current_Tag[P.STRINGLENGTH+1+P.PAIR_LENGTH_RIGHT]!='\n' && Current_Tag[P.STRINGLENGTH+1+P.PAIR_LENGTH_RIGHT]!='\r' && Current_Tag[P.STRINGLENGTH+1+P.PAIR_LENGTH_RIGHT]!=0;P.PAIR_LENGTH_RIGHT++);
			fgets(Quality,MAXTAG,Input_File);//plus
			if (Quality[0]=='>') P.FILETYPE=FA;else P.FILETYPE=FQ;
			if (P.FILETYPE == FQ && Quality[0] != '+' && Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
		}
		else
		{
			P.NORMAL_TAGS=TRUE;
			fgets(Quality,MAXTAG,Input_File);//plus
			if (Quality[0]=='>') P.FILETYPE=FA;else P.FILETYPE=FQ;
			if (P.FILETYPE == FQ && Quality[0] != '+' && Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
			fgets(Quality,MAXTAG,Input_File);//phred
		}

		fseek(Input_File, 0L, SEEK_END);
		P.File_Size = ftello64(Input_File);
		P.Org_File=Input_File;

		fseek(Input_File,0,SEEK_SET);//go top
		//if(!P.NORMAL_TAGS) fseek(Mate_File,0,SEEK_SET);//go top
	}
	READLEN=P.TAG_COPY_LEN;
	ORG_STRINGLENGTH=READLEN;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_Tag 
 *  Description:  Read next tag and store details in Head.. 
 * =====================================================================================
 */

char Read_Tag(READ & Head,READ & Tail,FILE* & Input_File,FILE* & Mate_File,FILETYPE & F,SAMREAD & SAM)
{
	char * Current_Tag;
	static int Random_Pointer=0;
	Current_Tag=Head.Tag;

	if(SAM.Enable) 
	{
		int Flag;
Read_Again:
		flockfile(Input_File);
		if (fgets(SAM.SAM_Line,5000,Input_File)!=0)// read a tag...
		{
			char NM[10];
			funlockfile(Input_File);
			if(*SAM.SAM_Line=='@')
				goto Read_Again;
			sscanf(SAM.SAM_Line,"%s %d %s %u %d %s %*s %*d %*d %s %s %s",Head.Description,&Flag,SAM.Chr,&SAM.Loc,&SAM.MapQ,SAM.Cigar,Head.Tag,Head.Quality,NM); 
			if(Flag!=4)
			{
				assert(*NM=='N');
				SAM.NM=atoi(NM+5);
			}
			else
				SAM.NM=INT_MAX;

			strcpy(Head.Tag_Copy,Head.Tag);
			Head.NCount=0;int j=0;
			for (unsigned i=0;i<=F.STRINGLENGTH-1;i++)
			{
				if (Current_Tag[i] == 'n' || Current_Tag[i]=='N')
				{
					Head.N[j++]=i;Head.NLocations[i]=TRUE;Head.NCount++;
					Current_Tag[i]=Random_Array[Random_Pointer++];Head.N[j++]=Current_Tag[i];
					if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
				}
				else Head.NLocations[i]=FALSE;
				if (Head.Tag_Copy[i]>='a') 
				{
					Head.Tag_Copy[i]-=('a'-'A'); 
				}
				Current_Tag[i]=Char_To_Code[Current_Tag[i]];
				Head.Complement[F.STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];
			}
			Current_Tag[F.STRINGLENGTH]='+';
			Head.Complement[F.STRINGLENGTH]='-';
			return TRUE;
		}
		else
		{
			funlockfile(Input_File);
			return FALSE;
		}
	}
	else if (fgets(Head.Description,MAXDES,Input_File)!=0)// read a tag...
	{
		fgets(Head.Tag,MAXDES,Input_File);//tag
		if (F.FILETYPE == FQ)
		{
			fgets(Head.Plus,MAXTAG,Input_File);//plus
			fgets(Head.Quality,MAXTAG,Input_File);//phred
		}
		funlockfile(Input_File);
		strcpy(Head.Tag_Copy,Head.Tag);
		Head.NCount=0;int j=0;
		for (unsigned i=0;i<=F.STRINGLENGTH-1;i++)
		{
			if (Current_Tag[i] == 'n' || Current_Tag[i]=='N')
			{
				Head.N[j++]=i;Head.NLocations[i]=TRUE;Head.NCount++;
				Current_Tag[i]=Random_Array[Random_Pointer++];Head.N[j++]=Current_Tag[i];
				if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
			}
			else Head.NLocations[i]=FALSE;
			if (Head.Tag_Copy[i]>='a') 
			{
				Head.Tag_Copy[i]-=('a'-'A'); 
			}
			Current_Tag[i]=Char_To_Code[Current_Tag[i]];
			Head.Complement[F.STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];
		}
		Current_Tag[F.STRINGLENGTH]='+';
		Head.Complement[F.STRINGLENGTH]='-';
	} 
	else 
	{
		funlockfile(Input_File);
		return FALSE;
	}
	if(F.NORMAL_TAGS) 
	{
		return TRUE;
	}
	

	Current_Tag=Tail.Tag;
	/*if (fgets(Mate_File,Tail.Description,MAXDES)!=0)// read a tag...
	{
		fgets(Mate_File,Tail.Tag,MAXDES);//tag
		if (F.FILETYPE == FQ)
		{
			fgets(Mate_File,Tail.Plus,MAXTAG);//plus
			fgets(Mate_File,Tail.Quality,MAXTAG);//phred
		}
		strcpy(Tail.Tag_Copy,Tail.Tag);
		Tail.NCount=0;int j=0;
		for (unsigned i=0;i<=F.STRINGLENGTH-1;i++)
		{
			if (Current_Tag[i] == 'n' || Current_Tag[i]=='N')
			{
				Tail.N[j++]=i;Tail.NLocations[i]=TRUE;Tail.NCount++;
				Current_Tag[i]=Random_Array[Random_Pointer++];Tail.N[j++]=Current_Tag[i];
				if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
			}
			else Tail.NLocations[i]=FALSE;
			Current_Tag[i]=Char_To_Code[Current_Tag[i]];
			Tail.Complement[F.STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];
		}
		Current_Tag[F.STRINGLENGTH]='+';
		Tail.Complement[F.STRINGLENGTH]='-';
		Current_Tag=Head.Tag;
		return TRUE;
	} else {printf("Read_Tag():Unpaired read...!\n");exit(0);};*/
}



void Open_Files(FILE* & Input_File,FILE* & Mate_File,Parameters P)
{
	Input_File=File_Open(P.PATTERNFILE,"r");//Load tags
	if(P.Patternfile_Count) Mate_File=File_Open(P.PATTERNFILE1,"r");//Load tags
}
