#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;
const int MAXTAG=200;
unsigned Genome_Size;

map <string,unsigned> Chromo_Info;
char Org_String[300];

void Load_Location(char* LOCATIONFILE, map <string,unsigned> & Chromo_Info);
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Bases
 *  Description:  read bases from packed file..
 * =====================================================================================
 */
void Get_Bases (unsigned Location,int StringLength,unsigned char* Original_Text, char* Org_String)
{
	Location--;
	for (int i=0;i<StringLength;i++)
	{
		unsigned char L= (unsigned char)(Original_Text[(Location+i)/4]<< (((Location+i) % 4) * 2)) >>6;
		Org_String[i]="ACGT"[L];
	}
}

int main(int argc,char* argv[])
{
	ifstream IN,BASE;
	stringstream SS;
	string input;
	unsigned char *Original_Text;
	unsigned count=0;
	try
	{
		if(argc !=4) throw("command line:\nlabel file length\n"); 
		/*IN.open(argv[1]);
		if(!IN.is_open()) throw ("cannot open BAM file...\n");
		BASE.open(argv[2],ios::binary);
		if(!BASE.is_open()) throw ("cannot open Base file...\n");
		else
		{
			BASE.seekg (0, ios::end);
			streampos Length=BASE.tellg();
			BASE.seekg (0, ios::beg);
			Original_Text=(unsigned char*) malloc(Length);
			BASE.read((char*) Original_Text,Length);
		}*/
		Load_Location(argv[1],Chromo_Info);
		int Loc=atoi(argv[3]);
		unsigned Flat_Location=Chromo_Info[argv[2]]+Loc;
		cout <<"Loc1 " <<argv[2]<<":"<<Loc<<endl;
		cout <<Genome_Size-Loc<<endl;
		cout <<Flat_Location<<endl;
		/*map <string,unsigned>::iterator I;
		while (getline(cin,input))
		{
			if(input.c_str()[0] != '@')
			{
				unsigned Flag,Loc;
				char Read_ID[MAXTAG];
				char Chr[MAXTAG];
				char Cigar[MAXTAG];
				char Seq[MAXTAG];
				sscanf(input.c_str(),"%s %u %s %u %*u %s %*s %*u %*u %s",Read_ID,&Flag,Chr,&Loc,Cigar,Seq);            
				if(!(Flag & 4))
				{
					int String_Len=strlen(Seq);
					int Cigar_Len=atoi(Cigar);
					//cout << input.c_str();
					if(Cigar_Len == String_Len) //no indels
					{
						unsigned Flat_Location=Chromo_Info[Chr]+Loc;
						Get_Bases (Flat_Location,String_Len,Original_Text,Org_String);
						//cout <<Read_ID<<"\t" <<Flag <<"\t" <<Chr<<"\t"  <<Loc<<"\t" <<Cigar<<"\t" <<Seq<<"\t"<<Org_String<<"\n";
						cout <<Read_ID<<"\t";
						int Mis=0;
						for (int i=0;i<String_Len;i++)
						{
							if(Seq[i] != Org_String[i])
							{
								SS <<i<<":"<<Seq[i]<<Org_String[i]<<"\t";
								Mis++;
							}
						}
						cout <<Mis<<"\t"<<SS.str()<<"\n";SS.str("");
					}
				}
			}
		}*/
	}
	catch(const char* err)
	{
		cout << err;
	}
}

//0       83      chr5    59992768        60      90M     =       59992054        -804    ACTTTCACTTGCCACAACCTCGTTGTGTGATCTTGCCATTCTAACTCTCAGATTCAATGAAAGCTGTATGTCTTTTAAAATTTATCAGTT   BBBBBBBBBBBBBBBBBBBBBBBBBBBBB_\_Yddce\ccddddd^baaddc^`deeeee\c\ffdffefcf`ffffeffee\eeee`eY      XT:A:U  NM:i:1  SM:i:37 AM:i:37      X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  MD:Z:4G8
/*int  Location_To_Genome(unsigned & Location,unsigned *Offsets,int Genome_Count)
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
}*/


void Load_Location(char* LOCATIONFILE, map <string,unsigned> & Chromo_Info)
{
	ifstream LOC;
	string Chr;
	unsigned Location,Tot_Location=0;

	LOC.open(LOCATIONFILE);
	if(!LOC.is_open()) throw ("Loc file not found..\n");
	//int Genome_Count=0;
	while (LOC>>Location && LOC>>Chr)
	{
		Tot_Location+=Location;
		Chromo_Info[Chr]=Tot_Location;
	}
	LOC>>Location;
	Genome_Size=Tot_Location+Location;
}
