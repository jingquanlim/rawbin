#include "Print.h"
extern char Char_To_CodeC[];
extern const int UNIQUE_SIGNAL;//there is one junction, it has a signal
extern const int UNIQUE_NOSIGNAL;//There is one junction, it does not have a signal..
extern const int NON_UNIQUE_SIGNAL;//There are many local junctions, but only one have a signal..
//extern Offset_Record Genome_Offsets[];
char const* Junc_Decode[]={"CAN","SEMI","SEMI","NON"};
extern bool PRINT_NON_CANON;
extern bool SOFTCLIP;
extern bool SAM_READER;
extern int READLEN;

void Reverse_Read(char* RDest,char* QDest,char* R,char* Q,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*RDest="ACGT"[Char_To_CodeC[R[i]]];RDest++;
		*QDest=Q[i];QDest++;
	}
	*QDest=0;
	*RDest=0;
}

char const* Get_Junc_Type(int Type)
{
	return Junc_Decode[3-Type];
	//return !(&& && 
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

char* Nullify_String(char* S)
{
	char* Strend;
	for(Strend=S;*Strend && *Strend!='\n';*Strend++);
	*Strend=0;
	return Strend;
}

void Print_Hits(READ & Head,Junction *Final_Juncs,ofstream & SAM,int firstSignal,unsigned Hit_ID,int Err,Offset_Record *Genome_Offsets,bool Multi_Hits,int READLENX,int Read_Skip)
{
	Nullify_String(Head.Description+1);
	Nullify_String(Head.Tag_Copy);
	Nullify_String(Head.Quality);
	assert(Hit_ID!=INT_MAX);//Hope there arent gazillion reads..
	char CIGAR[500];
	
	int i=firstSignal;
	int Off_Des=1;if(SAM_READER) Off_Des=0;
	if(Final_Juncs[i].q)
	{
		OP Pair;

		Ann_Info A;
		char* CIGAR_ptr=CIGAR;
		int Last=0;int Label=Final_Juncs[i].Label;int Bdr_Count=0;
		if(Read_Skip>0)
		{
			if(SOFTCLIP)
				CIGAR_ptr+=sprintf(CIGAR_ptr,"%dS",Read_Skip);
			else
			{
			Final_Juncs[i].r+=Read_Skip;
			READLENX+=Read_Skip;
			}
		}
		for(int j=0;j<Final_Juncs[i].Junc_Count;j++)
		{
			Location_To_Genome(Final_Juncs[i+j].p,A);//TODO:one lookup should be enough..
			Location_To_Genome(Final_Juncs[i+j].q,A);

			Pair.x=Final_Juncs[i+j].p;
			Pair.y=Final_Juncs[i+j].q;
			if(Final_Juncs[i+j].ID==INT_MAX-1) Final_Juncs[i].ID=A.ID;//Junction at extreme..
			assert(Genome_Offsets[Final_Juncs[i+j].ID].Offset !=INT_MAX && Final_Juncs[i+j].ID==Final_Juncs[i].ID );
			if(!Err && !Multi_Hits) Genome_Offsets[Final_Juncs[i].ID].Junc_Hash->Insert(Pair,Final_Juncs[i+j],Hit_ID? false:true);

			assert(Label==Final_Juncs[i].Label);
			CIGAR_ptr+=sprintf(CIGAR_ptr,"%dM%dN",Final_Juncs[i+j].r,Final_Juncs[i+j].q-Final_Juncs[i+j].p+1);
			Last+=Final_Juncs[i+j].r;
		}
		if(!SOFTCLIP)
		{
			if(Read_Skip<0)
				READLENX-=Read_Skip;
		}
		CIGAR_ptr+=sprintf(CIGAR_ptr,"%dM",READLENX-Last);
		if(SOFTCLIP)
		{
			if(Read_Skip<0)
				sprintf(CIGAR_ptr,"%dS",-Read_Skip);
		}

		char *R=Head.Tag_Copy,*Q=Head.Quality;char RQ[MAXTAG],RR[MAXTAG];
		if(!Final_Juncs[i].Sign)
		{
			Reverse_Read(RR,RQ,Head.Tag_Copy,Head.Quality,READLEN);
			Q=Head.Quality;R=Head.Tag_Copy;
		}

		SAM << Head.Description+Off_Des <<"\t"
			<< ((Final_Juncs[i].Sign) ? 0:16) <<"\t"  
			<< A.Name <<"\t" 
			<< Final_Juncs[i].p-Final_Juncs[i].r+1 << "\t" 
			<< (int)((Multi_Hits)? 0:60) <<"\t" 
			<<CIGAR<<"\t"
			<< "*\t0\t0\t" 
			<< R << "\t" 
			<< Q << "\t"
			<< "XH:i:" << Final_Juncs[i].Mismatches <<"\t" << "XJ:i:" << Final_Juncs[i].score 
			<<"\tXS:A:"<<(char) Canonical_Sign(Final_Juncs[i].signal)<<"\tER:i:"<<Err<<endl;
	}
	else//Full read map..
	{
		assert(Final_Juncs[i].r==0);
		Ann_Info A;
		Location_To_Genome(Final_Juncs[i].p,A);
		char *R=Head.Tag_Copy,*Q=Head.Quality;char RQ[MAXTAG],RR[MAXTAG];
		if(!Final_Juncs[i].p || (Final_Juncs[i].p+READLEN > A.Size))//check for a Boundary Hit..
			return;
		if(!Final_Juncs[i].Sign)
		{
			Reverse_Read(RR,RQ,Head.Tag_Copy,Head.Quality,READLEN);
			Q=RQ;R=RR;
		}
/////////////////////////////////////////////////////////
		if(!SOFTCLIP)
		{
			if(Read_Skip>0)
			{
				Final_Juncs[i].p-=Read_Skip;
				READLENX+=Read_Skip;
				if(Final_Juncs[i].p ==0) 
					return;
			}
			else if(Read_Skip<0)
			{
				READLENX-=Read_Skip;
			}
		}
/////////////////////////////////////////////////////////
		SAM << Head.Description+Off_Des <<"\t"
			<< ((Final_Juncs[i].Sign) ? 0:16) <<"\t"  /*Flag*/
			<< A.Name <<"\t" 
			<< Final_Juncs[i].p << "\t" 
			<< (int)((Multi_Hits)? 0:60) <<"\t";

		if(SOFTCLIP)
		{
			if(Read_Skip>0)
				SAM << Read_Skip<<"S";
			SAM << READLENX << "M";  
			if(Read_Skip<0)
				SAM << -Read_Skip<<"S";
		}
		else
			SAM << READLENX << "M";  
		SAM	<< "\t*\t0\t0\t" 
			<< R << "\t" 
			<< Q << "\t"
			<< "XH:i:" << Final_Juncs[i].Mismatches <<"\t" <<"XS:A:." <<"\tER:i:"<<Err<<endl;
	}
}

void Print_SAM_Header(std::map <unsigned, Ann_Info> Annotations,int argc,char* argv[],char* Input_File)
{
	ofstream SAMHEAD;
	SAMHEAD.open("0header.sam");
	if(SAMHEAD.is_open())
	{
		std::map <unsigned, Ann_Info> ::iterator S,E;
		S=Annotations.begin();E=Annotations.end();
		unsigned CSize=0;
		while (S!=E)
		{
			Ann_Info T=S->second;
			if(T.Size > CSize) CSize=T.Size;
			SAMHEAD<<"@SQ\tSN:"<<T.Name<<"\tLN:"<<T.Size<<endl;
			S++;
		}
		char Current_Dir[1000];
		if (!getcwd(Current_Dir,990))
		{
			sprintf (Current_Dir,"%d",rand());
		}
		SAMHEAD<<"@RG\tID:"<<Current_Dir<<"\tSM:"<< Input_File <<endl;
		SAMHEAD<<"@PG\tID:PEnGuin\tCL:";
		for(int i = 0; i < argc; i++) SAMHEAD << " "<<argv[i];
		SAMHEAD << endl;
	}
	else
	{
		cout << "Print_SAM_Header():Header cannot be opened..\n";
		exit(100);
	}
}


void Print_Junctions(char* Junction_File,Offset_Record *Genome_Offsets)
{
	static char NOHEADER=TRUE;
	ofstream JUNCFILE;
	JUNCFILE.open(Junction_File);
	if(!JUNCFILE.is_open())
	{
		cout << "Print_Junctions(): Junction file cannot be opened ..\n";
		exit(100);
	}

	for(int i=0;Genome_Offsets[i].Offset !=INT_MAX;i++)
	{
		Hash *Junctions=Genome_Offsets[i].Junc_Hash;
		if (NOHEADER) 
		{
			JUNCFILE <<"track name=\"Rawbin Exons\"\n";
			if(JUNCFILE.fail()) { cout << "Print_Junctions(): Error writing to junction file\n";exit(100);} 
			NOHEADER=FALSE;
		}
		char* Chromosome=Genome_Offsets[i].Genome;
		cout << "Dumping junctions " << Chromosome <<endl;	
		static int j=0;
		JStat JStat;
		OP JPair;
		char Junc_Not_Empty = Junctions->Init_Iterate(JPair,JStat);

		while (Junc_Not_Empty)
		{
			//if(JStat.Junc_Type)
			if(JStat.Unique && (JStat.Junc_Type || PRINT_NON_CANON))
			{

				OP H,T;

				char Strand_Sign = '-';//(JStat.Junc_Type >= MINUS_JUNC) ? '-' : '+';
				int Left_Cover=100;//Coverage[JPair.x];
				int Right_Cover=100;//Coverage[JPair.y];
				int Left_Drop=10;//Coverage[JPair.x]-Coverage[JPair.x+1];
				int Right_Drop=100;//Coverage[JPair.y]-Coverage[JPair.y-1];

				//H.x=JPair.x-2;H.y=JPair.x+1;
				//T.x=JPair.y;T.y=JPair.y+2;
				JPair.x=JPair.x-1;
				JPair.y=JPair.y+1;
				H.x=JPair.x-2;H.y=JPair.x+1;
				T.x=JPair.y;T.y=JPair.y+2;
				if(H.x>=T.x || T.x >= T.y || H.x >= H.y)
				{
					printf ("Print_BED():Junction Parse error..\nH.x:H.y=[%u,%u]       T.x:T.y=[%u,%u]\n",H.x,H.y,T.x,T.y);exit(0);
				}
				char const *Junction_Type=Get_Junc_Type(JStat.Junc_Type);
				JUNCFILE
					<< Chromosome <<"\t"
					<< H.x <<"\t" << T.y+1 <<"\t"/*Chrom start and end...*/
					<< "JUNC"
					<< j++ <<"_" /*Junction name..*/
					<< Junction_Type <<"_"
					<< Left_Cover <<"_"
					<< Right_Cover <<"_"
					<< Left_Drop <<"_"
					<< Right_Drop <<"\t"
					<< JStat.Count <<"\t"/*score to color...*/
					<< Strand_Sign <<"\t"

					<< H.x <<"\t" /*thickStart and end*/
					<< T.y <<"\t"
					<<"255,0,0\t2\t"

					<< H.y-H.x <<","
					<< T.y-T.x+1 <<"\t"/*block sizes..*/
					<< 0 << ","
					<< T.x-H.x//block starts...
					<< endl;
				if(JUNCFILE.fail()) { cout << "Print_Junctions(): Error writing to junction file\n";exit(100);} 

			}
			Junc_Not_Empty=Junctions->Iterate(JPair,JStat);
		}
	}
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

/*void Print_BED(Hash & Junctions,Island* Island_List, COVTYPE* Coverage,char* chromosome,char Strand)
{
	static char NOHEADER=TRUE;
	printf("\t+Writing Files..\n");

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
/*			if(Pass_Junc)
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
}*/
