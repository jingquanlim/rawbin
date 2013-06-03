#include "init.h"
extern int INDEX_RESOLUTION;
extern int EXONGAP;
extern int READLEN;
extern int RQFACTOR;
extern unsigned CONVERSION_FACTOR;
extern char Char_To_CodeC[256];
extern char Char_To_Code[256];
extern unsigned RQ_Hits;
extern SA* SA_Index;
extern char COMPRESS;
extern int *SA_Blocks;
extern std::map <unsigned, Ann_Info> Annotations;
extern unsigned Location_Array[80];
extern pthread_mutex_t OpenGenomeFileslock;
extern bool SAM_READER;
extern int ORG_STRINGLENGTH; 
void Init(BWT *revfmi,unsigned & SOURCELENGTH,FILE* & Input_File,FILE* & Mate_File,FILETYPE & File_Info,Parameters & CL,Index_Info & Genome_Files,int & INIT_MIS_SCAN)
{
	EXONGAP=CL.EXONGAP;
	SOURCELENGTH = revfmi->textLength;
	CONVERSION_FACTOR=revfmi->textLength-RQFACTOR;//+1;
	Char_To_Code['N']=0;Char_To_Code['n']=0;
	Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;
	Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;
	Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;
	Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;
	Char_To_CodeC['A']=3;Char_To_CodeC['C']=2;Char_To_CodeC['G']=1;Char_To_CodeC['T']=0;
	Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;
	Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt

	Open_Files(Input_File,Mate_File,CL);
	if(SAM_READER)
	{
		char SAM_Line[5000];
		int Flag;READ Head;SAMREAD SAM;
Read_Again:
		if (fgets(SAM_Line,5000,Input_File)!=0)// read a tag...
		{
			if(*SAM_Line=='@')
				goto Read_Again;
			sscanf(SAM_Line,"%s %d %s %u %d %s %s %*d %*d %s ",Head.Description,&Flag,SAM.Chr,&SAM.Loc,&SAM.MapQ,SAM.Cigar,Head.Quality,Head.Tag_Copy); 
		}
		else
		{
			cout << "Init():Error reading file..\n";exit (100);
		}
		File_Info.STRINGLENGTH=strlen(Head.Tag_Copy);
		if (Head.Quality[0]=='*')
		{
		       	File_Info.FILETYPE=FA;
		}
		else 
		{
			assert(strlen(Head.Quality)==File_Info.STRINGLENGTH);
			File_Info.FILETYPE=FQ;
		}

		fseek(Input_File, 0L, SEEK_END);
		File_Info.File_Size = ftello64(Input_File);
		File_Info.Org_File=Input_File;
		fseek(Input_File,0,SEEK_SET);//go top

		READLEN=File_Info.TAG_COPY_LEN=File_Info.STRINGLENGTH;
		ORG_STRINGLENGTH=READLEN;
	}
	else
		Detect_Input(File_Info,Input_File,Mate_File);
	FILE* Inf_File;
	if((Inf_File=File_Exist_Open(Genome_Files.INFOFILE)))
	{
		char String[40];
		int Num;
		fscanf(Inf_File,"%s%d",String,&Num);
		fscanf(Inf_File,"%s%d",String,&INDEX_RESOLUTION);
	}
	else
	{
		INDEX_RESOLUTION=30000;
	}
	
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

	//fclose(Inf_File);
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

void Load_All_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,RANGEINDEX & Range_Index)
{
	if(mkdir("Raw_Out",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) && errno != EEXIST) {fprintf(stderr,"Cannot create Temp directory..\n");exit(-1);};
	fprintf (stderr,"Loading index %s\n",Genome_Files.BWTFILE);
	fwfmi= Load_Indexes(Genome_Files.BWTFILE,Genome_Files.OCCFILE,Genome_Files.SAFILE,mmPool);
	fprintf (stderr,"Loading index %s\n",Genome_Files.REVBWTINDEX);
	revfmi= Load_Indexes(Genome_Files.REVBWTINDEX,Genome_Files.REVOCCFILE,Genome_Files.REVSAFILE,mmPool);
	fwfmi->saInterval=revfmi->saInterval;
	fprintf (stderr,"Loading Location %s\n",Genome_Files.LOCATIONFILE);
	int Genome_Count= Load_Location(Genome_Files.LOCATIONFILE,Annotations,Location_Array);
	fprintf (stderr,"Loading packed genome %s\n",Genome_Files.BINFILE);
	loadPac(Genome_Files.BINFILE);
	fprintf (stderr,"Loading index %s\n",Genome_Files.INDFILE);
	Load_Range_Index(Genome_Files.INDFILE,Genome_Files.BLKFILE,Range_Index);
	SA_Index=Range_Index.SA_Index;
	SA_Blocks=Range_Index.SA_Blocks;
	COMPRESS=Range_Index.COMPRESS;
	RQ_Hits=Range_Index.Hits;
	fprintf(stderr,"Done...\n");
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Open_Genome_Files
 *  Description:  Creates an array of existing chromosomes and associates auxilliary files
 *  		  and data structures with it.
 * =====================================================================================
 */
int Open_Genome_Files(char* LOCATIONFILE,Offset_Record* Genome_Offsets,unsigned Offsets[])
{
	pthread_mutex_lock(&OpenGenomeFileslock);
	FILE* Location_File=File_Open(LOCATIONFILE,"r");
	int Genome_Count=0;

	/*if(mkdir("Raw_Out",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) && errno != EEXIST) 
	{
		printf("Open_Genome_Files():Cannot create Temp directory..\n");
		exit(-1);
	}*/

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
	fclose(Location_File);
	for ( int i=1;i<Genome_Count;i++)
	{
		Offsets[i]=Offsets[i-1]+Genome_Offsets[i].Offset;
		//if(MAPMODE) 
		{
			Genome_Offsets[i-1].Junc_Hash=new Hash;
			/*string S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_File=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_FileM=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].Unmapped=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].UnmappedM=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedX=File_Open(S.c_str(),"w+b");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedXM=File_Open(S.c_str(),"w+b");*/
			
		}
		/*else
		{
			string S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_File=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S;Genome_Offsets[i-1].Out_FileM=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].Unmapped=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".U";Genome_Offsets[i-1].UnmappedM=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].Genome;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedX=File_Open(S.c_str(),"r");
			S=Genome_Offsets[i-1].GenomeM;S="Raw_Out/"+S+".UX";Genome_Offsets[i-1].UnmappedXM=File_Open(S.c_str(),"r");
		}*/
		/*if (USEREFGENE) 
		{
			string S=Genome_String+Genome_Offsets[i-1].Genome+".ref";
			Genome_Offsets[i-1].Ref_File=File_Open(S.c_str(),"r");
		}*/
	}
	Genome_Offsets[Genome_Count-1].Offset=INT_MAX;
	pthread_mutex_unlock(&OpenGenomeFileslock);
	return Genome_Count-1;
}


void Load_FM_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool)
{
	fprintf (stderr,"Loading index %s\n",Genome_Files.BWTFILE);
	fwfmi= Load_Indexes(Genome_Files.BWTFILE,Genome_Files.OCCFILE,Genome_Files.SAFILE,mmPool);
	fprintf (stderr,"Loading index %s\n",Genome_Files.REVBWTINDEX);
	revfmi= Load_Indexes(Genome_Files.REVBWTINDEX,Genome_Files.REVOCCFILE,Genome_Files.REVSAFILE,mmPool);
	fwfmi->saInterval=revfmi->saInterval;
	fprintf (stderr,"Loading Location %s\n",Genome_Files.LOCATIONFILE);
	int Genome_Count= Load_Location(Genome_Files.LOCATIONFILE,Annotations,Location_Array);
	/*fprintf (stderr,"Loading packed genome %s\n",Genome_Files.BINFILE);
	loadPac(Genome_Files.BINFILE);*/
	fprintf (stderr,"Loading index %s\n",Genome_Files.INDFILE);
	fprintf(stderr,"Done...\n");
}
