#ifndef __HASH_GUARD__
#define __HASH_GUARD__
#include <map>
#define ENT_LIM 50
using namespace std;

struct Junction {
	unsigned p;
	unsigned q;
	unsigned r;
	char *Chrom;
	int Mismatches;
	char signal[5];
	int isCanonical();
	int Label;
	int ID;
	char Sign;
	float score;
	int Type;
	int R,L;
	int Junc_Count;
};


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
	bool Unique;
	char L_Anchor[ENT_LIM];
	char R_Anchor[ENT_LIM];
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

typedef map <OP,JStat>::iterator map_it; 
class Hash
{

	public:
	map <OP,JStat,OP_Cmp> Junctions;
	map_it Junc_I,JJ;
	map_it Last;

	void Insert (OP Location,Junction & J,bool Unique);
	void Delete (OP Location);
	bool Init_Iterate(OP & Location,JStat & Data);
	bool Iterate(OP & Location,JStat & Data);
	map_it Begin();


};
#endif
