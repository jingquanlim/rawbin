#include "extend.h"
#include "limits.h"
#include "file.h"
#include "batlib.h"
#include "common.h"
#include <assert.h>
using namespace std;

extern const int MIS_DENSITY;
extern const bool DEBUG;
extern const int MINX;//Minimum extension..
extern float Donor_Prob[16][64][2];
extern float Acc_Prob[16][64][2];

bool Junction::isCanonical(){
	return !(strcmp(signal, "GTAG") && strcmp(signal, "GCAG") && strcmp(signal, "ATAC") && strcmp(signal, "CTAC") && strcmp(signal, "CTGC") && strcmp(signal, "GTAT"));
}

/*int main(int argc, char* argv[]){
	loadPac("test.fasta.pac");
	//printf("%s\n",Original_Text);
	char* R = "GCATCGATCAGCATGCATCGATGGCATGCGTGGGAGAATGTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
	Junction* junctions = extend(R, 20, 100, 560,680);
	
	delete [] junctions;
	return 0;
}*/

void loadPac(char* filename){
	FILE* Original_File=File_Open(filename,"rb");
	Original_Text=(unsigned char*) malloc(Get_File_Size(Original_File));
	fread(Original_Text,Get_File_Size(Original_File),1,Original_File);
}

void Get_Bases_ASCII (unsigned Location,int StringLength,char* Org_String)
{
	for (int i=0;i<StringLength;i++)
	{
		unsigned char L= (unsigned char)(Original_Text[(Location+i)/4]<< (((Location+i) % 4) * 2)) >>6;
		Org_String[i]="ACGT"[L];
	}
	Org_String[StringLength]=0;
}

void Get_Bases(unsigned Location,int StringLength,char* Org_String)
{
	for (int i=0;i<StringLength;i++)
	{
		unsigned char L= (unsigned char)(Original_Text[(Location+i)/4]<< (((Location+i) % 4) * 2)) >>6;
		Org_String[i]=L;
	}
}

float signalScore(char* signal) {
	//printf("%s\n",signal);
	if(strcmp(signal,"GTAG") == 0 || strcmp(signal,"CTAC")==0)
		return 4;
	else if(strcmp(signal,"GCAG")==0 || strcmp(signal,"CTGC")==0)
		return 2;
	else if(strcmp(signal,"ATAC")==0 || strcmp(signal,"GTAT")==0)
		return 1;
	else
		return 0;
}

inline void Convert_Reverse(char* Read_bin,char * RC_bin,int StringLength)
{
	for (unsigned i=0;i<=StringLength-1;i++)
	{
		RC_bin[StringLength-1-i]="ACGT"[3-Read_bin[i]];
	}
}

inline void b2c(char* bin, char* read, int length) {
	for(int i =0;i<length; i++) {
		read[i] = "ACGT"[bin[i]]; 
	}
}

float getScore(unsigned p, unsigned q, unsigned x, int size, int i, int misL, int misR, char sign) {
	float score = 0-misL-misR;
	
	//closed interval, 0-based
	unsigned juncStart = p+x+i+1;
	unsigned juncEnd = q-size+i-1;
	
	//if(juncEnd+1 == juncStart) {
	//	score += 5;
	//	return score;
	//}
	char signal[5];

	//Get the score contributed by signal.
	Get_Bases_ASCII(juncStart,2,signal);
	Get_Bases_ASCII(juncEnd-1,2,signal+2);
	score += signalScore(signal);
	
	char donor[8],acc[8];
	//Get the score of the surrounding sequence.
	if(sign == '+') {
		Get_Bases(juncStart-3, 8, donor);
		Get_Bases(juncEnd-4, 8, acc);
	} else {
		char donor_temp[8],acc_temp[8];
		Get_Bases(juncStart-3, 8, acc_temp);
		Get_Bases(juncEnd-4, 8, donor_temp);
		Convert_Reverse(acc_temp, acc,8);
		Convert_Reverse(donor_temp, donor,8);
	}
	float donor_ex_prob = Donor_Prob[4*donor[3]+donor[4]][16*donor[0]+4*donor[1]+donor[2]][0];
	float donor_in_prob = Donor_Prob[4*donor[3]+donor[4]][16*donor[5]+4*donor[6]+donor[7]][1];
	float acc_ex_prob = Acc_Prob[4*acc[3]+acc[4]][16*acc[0]+4*acc[1]+acc[2]][0];
	float acc_in_prob = Acc_Prob[4*acc[3]+acc[4]][16*acc[5]+4*acc[6]+acc[7]][1];
	score += std::min(donor_ex_prob,donor_in_prob)*std::min(acc_ex_prob,acc_in_prob);
	
	char donorStr[9],accStr[9];
	donorStr[8] = 0;
	accStr[8] = 0;
	b2c(donor,donorStr,8);
	b2c(acc,accStr,8);
	Ann_Info A, A1;
	Location_To_Genome(juncStart,A);
	Location_To_Genome(juncEnd,A1);
	Location_To_Genome(p, A1);
	Location_To_Genome(q, A1);
	printf("%s\t%u\t%u\t%d\t%d\t%u\t%u\t%s\t%s\t%f\t%d\t%d\t%c\t%s\n",A.Name,p,q,size,i,juncStart, juncEnd,donorStr,accStr,score,misL,misR,sign,signal);

	return score;
}

//Assume the coordinates are 0-based.
//All intervals used are closed interval.
Junction* extend(char* R, unsigned x, unsigned y, unsigned p, unsigned q, char sign){
	assert(p+x < q);
	
	//To allow some over extension into the anchors
	x -= 5;
	y += 5;
	q += 5;


	unsigned size = y-x-1;
	//printf("size: %u\n",size);
	int misL[size+1], misR[size+1];
	char basesL[size+1], basesR[size+1];
	basesL[size] = 0;
	basesR[size] = 0;
	//char str1[602];
	//str1[600]=0;
	//
	
	//printf("x: %u, y: %u, p: %u, q: %u\n",x,y,p,q);

	Get_Bases_ASCII(p+x+1, size, basesL);
	Get_Bases_ASCII(q-size, size, basesR);
	//Get_Bases_ASCII(0, 600, str1);
	//printf("size: %d, %d, %s\n",size,q-size, str1);
	//printf("%s\n",R);
	//printf("%c\t%d\t%s\n",basesL[0],strlen(basesL),basesL);
	//printf("%c\t%d\t%s\n",basesR[0],strlen(basesR),basesR);
	int countL=0, countR=0;
	misL[0] = misR[size] = 0;
	for(int i=1; i<size+1; i++)
	{
		if(R[x+i] != basesL[i-1])
			++countL;
		if(R[x+size+1-i] != basesR[size-i])
			++countR;
		//printf("%c\t%c\n",R[x+i], basesL[i-1]);
		misL[i] = countL;
		misR[size-i] = countR;
		//printf("misL: %d, %d\n",i,misL[i]);
		//printf("misR: %d, %d\n",i,misR[size-i]);
	}

	//Find the partitions that give the minimum mismatches.
	//Partitions store the size of the portion of the left extension.
	int min= misL[0] + misR[0];
	//float max = getScore(p,q,x,size,0,misL[0],misR[0],sign);
	int parCount = 0;
	int partitions[size+1];
	partitions[0] = 0;
	for(int i=1; i<size+1; i++)
	{
		//printf("%d\n",min);
		int temp = misL[i]+misR[i];
		//float temp = getScore(p,q,x,size,i,misL[i],misR[i],sign);
		if(DEBUG)
		{
			if(temp<=MIS_DENSITY)
			{
				partitions[parCount] = i;
				++parCount;
			}
			min=std::min(temp,min);
		}
		else
		{
			
			if(temp == min)//(fabs(temp-max) < 0.00000001)//(temp == min)
			{
				partitions[parCount] = i;
				++parCount;
			}
			else if(temp < min) //(temp>max)
			{
				parCount = 1;
				min = temp;
				//max = temp;
				partitions[0] = i;
			}
		}
	}
	//printf("Parcount: %d\t%f\n", parCount, max);
	
	//Find the junctions based on the partitions obtained
	//assert(min>=0 && min<=MINX+1);
	Junction *junctions = new Junction[parCount+1];
	for(int i = 0; i<parCount; i++)
	{
		junctions[i].p = p + x + partitions[i] + 1;
		junctions[i].q = q - size + partitions[i] - 1;
		junctions[i].r = x + partitions[i] + 1;
		junctions[i].score = getScore(p,q,x,size,partitions[i],misL[partitions[i]],misR[partitions[i]],sign);
		//junctions[i].Mismatches = min;
		junctions[i].Mismatches = misL[partitions[i]] + misR[partitions[i]];
		char signal[5];
		signal[4] = 0;
		Get_Bases_ASCII(junctions[i].p,2,signal);
		Get_Bases_ASCII(junctions[i].q-1,2,signal+2);
		//strcpy(donor,"GTAG");
		strcpy(junctions[i].signal, signal);
		//printf("%d\t%d\t%d\t%s\n",junctions[i].p,junctions[i].q,junctions[i].r,junctions[i].signal);
		//if(junctions[i].isCanonical())
		//	printf("canonical\n");
		//else
		//	printf("non-canonical\n");
	}
	junctions[parCount].p = UINT_MAX;

	return junctions;
}
