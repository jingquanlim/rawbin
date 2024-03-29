#ifndef __INIT_GUARD__
#define __INIT_GUARD__
#include <stdio.h>
#include <string>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include "bfix.h"
#include <getopt.h>
#include "zlib.h"
#include "const.h"
#include "assert.h"
#include "Hash.h"
#include "rqindex.h"
#include <fstream>
#include <iostream>
#include "batlib.h"
#include "Cmdline.h"
#include "Indexes.h"
#include "file.h"
#include "extend.h"
#include <pthread.h>

void Init(BWT *revfmi,unsigned & SOURCELENGTH,FILE* & Input_File,FILE* & Mate_File,FILETYPE & File_Info,Parameters & CL,Index_Info & Genome_Files,int & INIT_MIS_SCAN);
bool  Progress_Bar(Parameters & CL,unsigned & Number_of_Tags,unsigned & Progress,unsigned & Tag_Count,FILETYPE & File_Info);
void Load_All_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,RANGEINDEX & Range_Index);
int Open_Genome_Files(char* LOCATIONFILE,Offset_Record* Genome_Offsets,unsigned Offsets[]);
void Load_FM_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool);

#endif
