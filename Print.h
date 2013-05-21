#ifndef __PRINT_GUARD__
#define __PRINT_GUARD__
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
void Print_Junctions(char* Junction_File,Offset_Record *Genome_Offsets);
void Open_Outputs(ofstream & SAM,string filename);
char* Nullify_String(char* S);
void Print_Hits(READ & Head,Junction *Final_Juncs,ofstream & SAM,int firstSignal,unsigned Hit_ID,int Err,Offset_Record *Genome_Offsets,bool Multi_Hits,int READLEN,int Read_Skip);
void Print_SAM_Header(std::map <unsigned, Ann_Info> Annotations,int argc,char* argv[],char* Input_File);
#endif
