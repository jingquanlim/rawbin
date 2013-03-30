
#ifndef ORDPAIR
#define ORDPAIR
#include <map>
using namespace std;
#define INIT 100
const int REFMARK=10000;

struct OP//Ordered pairs..
{
	unsigned x;
	unsigned y;
};

struct JStat//Junction stats..
{
	int Count;
	int Gap;
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

class Hash
{

	public:
	map <OP,JStat,OP_Cmp> Junctions;
	map <OP,JStat,OP_Cmp>::iterator Junc_I;
	map <OP,JStat,OP_Cmp>::iterator Last;

	void Insert(OP Location,int Gap);//insert ordered pair...
	void Insert (OP Location,int Paring,int Gap);
	void Insert_Ref (OP Location,int Paring,int Gap);
	void Delete (OP Location);
	bool Init_Iterate(OP & Location,JStat & Data);
	bool Iterate(OP & Location,JStat & Data);
	bool Init_IterateX(OP & Location,JStat & Data);
	bool IterateX(OP & Location,JStat & Data);

};
#endif
