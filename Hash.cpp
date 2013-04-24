#include "Hash.h"
#include "limits.h"
#include "assert.h"
#include <iostream>
extern const int CANONICAL_SCORE;

map_it Hash::Begin ()
{
	return Junctions.begin();
}

void Hash::Insert (OP Location,Junction & J,bool Unique)
{
	JStat T;
	int Paring=J.Type;
	assert(Paring >=0 && Paring <=CANONICAL_SCORE); 
	assert(Location.x <= Location.y);
	Junc_I=Junctions.find(Location);
	if (Junc_I == Junctions.end()) //New entry..
	{
		T.Count=1;T.Junc_Type=Paring;
		/*for(int i=0;i<ENT_LIM-15;i++)
		{
			T.L_Anchor[i]=0;
			T.R_Anchor[i]=0;
		}*/
		if (Unique)
			T.Unique=true;
		else
			T.Unique=false;
		Junctions[Location]=T;
	}
	else 
	{
		//assert(J.L>=0 && J.R>=0);
		assert((Junc_I->second).Count>0 && (Junc_I->second).Junc_Type==Paring);
		(Junc_I->second).Count++;
		/*if(J.L>=ENT_LIM)
		{
			if ((Junc_I->second).L_Anchor[ENT_LIM-1]!=CHAR_MAX) (Junc_I->second).L_Anchor[ENT_LIM-1]++;
		}
		else if ((Junc_I->second).L_Anchor[J.L]!=CHAR_MAX) (Junc_I->second).L_Anchor[ENT_LIM-1]++;

		if(J.L>=ENT_LIM)
		{
			if ((Junc_I->second).R_Anchor[ENT_LIM-1]!=CHAR_MAX) (Junc_I->second).R_Anchor[ENT_LIM-1]++;
		}
		else if ((Junc_I->second).R_Anchor[J.L]!=CHAR_MAX) (Junc_I->second).R_Anchor[ENT_LIM-1]++;*/
		if (Unique)
			(Junc_I->second).Unique=true;
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
