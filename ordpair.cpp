#include "ordpair.h"
void Hash::Insert (OP Location,int Gap)
{
	Junc_I=Junctions.find(Location);
	if (Junc_I == Junctions.end()) //New entry..
	{
		(Junctions[Location]).Count=1;
		(Junctions[Location]).Gap=Gap;
		(Junctions[Location]).Junc_Type=INIT;
	}
	else 
	{
		(Junc_I->second).Count++;//=Data;
		(Junc_I->second).Gap=Gap;
	}
}

void Hash::Insert (OP Location,int Paring,int Gap)
{
	Junc_I=Junctions.find(Location);
	if (Junc_I == Junctions.end()) //New entry..
	{
		(Junctions[Location]).Count=1;
		(Junctions[Location]).Gap=Gap;
		(Junctions[Location]).Junc_Type=Paring;
	}
	else 
	{
		(Junc_I->second).Count++;//=Data;
		(Junc_I->second).Gap=Gap;
		(Junc_I->second).Junc_Type=Paring;
	}
}

void Hash::Insert_Ref (OP Location,int Paring,int Gap)
{
	Junc_I=Junctions.find(Location);
	if (Junc_I == Junctions.end()) //New entry..
	{
		(Junctions[Location]).Count=REFMARK;
		(Junctions[Location]).Gap=Gap;
		(Junctions[Location]).Junc_Type=Paring;
	}
	else 
	{
		(Junc_I->second).Count++;//=Data;
		(Junc_I->second).Gap=Gap;
		(Junc_I->second).Junc_Type=Paring;
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

bool Hash::Init_IterateX(OP & Location,JStat & Data)
{
	Junc_I=Junctions.begin();
	while(!((Junc_I->second).Junc_Type)) if (Junc_I == Junctions.end()) break; else Junc_I++;//skip deleted / no juction cases..
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

bool Hash::IterateX(OP & Location,JStat & Data)
{
	if(Junc_I == Junctions.end()) return false;
	while(!((Junc_I->second).Junc_Type)) if (Junc_I == Junctions.end()) return false; else Junc_I++;//skip deleted / no juction cases..
	Last=Junc_I;//save last pos...
	Location = Junc_I->first;
	Data = Junc_I->second;
	Junc_I++;
	return true;
}

