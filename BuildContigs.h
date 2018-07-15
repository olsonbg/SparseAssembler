#ifndef __BUILD_CONTIGS_H
#define __BUILD_CONTIGS_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include "BasicDataStructure.h"
#include "GraphConstruction.h"

//mark the branches in the sparse kmer graph
void MarkBranches(hashtable *ht);

void MarkBranches2(hashtable2 *ht);

void MarkBranches3(hashtable3 *ht);

void MarkBranches4(hashtable4 *ht);


void MarkBranches0(hashtable0 *ht);

// produce the contigs, single end assembly complete
void build_contigs(struct hashtable *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips);

void build_contigs2(struct hashtable2 *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips);


void build_contigs3(struct hashtable3 *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips);

void build_contigs4(struct hashtable4 *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips);

void build_contigs0(struct hashtable0 *ht, key_table *key_table, int K_size, int gap,string Contig_Filename,bool ScreenOffTips);

#endif

