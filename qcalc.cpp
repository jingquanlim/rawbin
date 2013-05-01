#include "qcalc.h"
extern float POW10[];
extern int POWLIMIT;
extern float QLIMIT_FLOAT;
extern float QLIMIT; 
extern int QUALITYCONVERSIONFACTOR;

float Pr(float Q)
{
	assert((int)Q>=0 && (int)Q<POWLIMIT-1);
	return(1-POW10[(int)Q]);
	//printf("table: %f\tlib: %f\n",POW10[(int)Q],1-pow(10,-Q/10));
	//return (1-pow(10,-Q/10));
}

float Pow10(float Q)
{
	static float Max=0;
	//if(Max<Q) {Max=Q;printf("Newmax: %f\n",Max);}
	assert((int)Q>=0);
	if((int)Q<POWLIMIT-1)
		return POW10[(int)Q];
	else
		return pow(10,-float(Q)/10);
}

void Build_Pow10()
{
	for(int Q=0;Q<POWLIMIT;Q++)
	{
		POW10[Q]=(pow(10,-float(Q)/10));
	}
}


float Penalty(char Q)
{
	float Q_Value=Q-QUALITYCONVERSIONFACTOR;//make quality to integer..
	//Q_Value=std::min(30.0f,Q_Value);
	assert(Q_Value>=0);// && Q_Value<=40);
	float Penalty= -10*log10((1-Pr(Q_Value))/3);
	Penalty=std::min(QLIMIT_FLOAT,Penalty);
	return Penalty;
}
