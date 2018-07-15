#ifndef __SCAFFOLDING_DATA_STRUCTURE_H
#define __SCAFFOLDING_DATA_STRUCTURE_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include "BasicDataStructure.h"

//The below is for scaffolding, can be omitted first %{

struct contig_t
{
	uint64_t contig_bits[100000];
	int contigLen;
};


struct adjacent_contig_info
{
	int32_t dist_sum;
	int8_t cov;
	string bridge;

};


struct scaffold_contig_info
{
	int32_t dist_sum;
	int8_t cov;

};

struct c_info
{

	uint8_t used:1,removed:1,unique:1,flip:1,marked:1,loc_unique:1,positive:1,rcomp:1;
};

//%} The above is for scaffolding, can be omitted first



bool it_cmp( const vector<int>::iterator &a,const vector<int>::iterator &b);

void Init_Contig(string &seq,struct contig_t & contig);

void ContigsRemapping(struct hashtable *ht,struct hashtable2 *ht2, int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt);

void ContigsRemapping3(struct hashtable3 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt);

void ContigsRemapping4(struct hashtable4 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt);

void ContigsRemapping0(struct hashtable0 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt);

void SuperContigsRemapping(struct hashtable *ht,struct hashtable2 *ht2, int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt);

void SuperContigsRemapping3(struct hashtable3 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt);

void SuperContigsRemapping4(struct hashtable4 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt);

void BuildContigAdjacency(hashtable *ht1, hashtable2 *ht2, struct contigs_info *contigs_info,int K_size, string ContigFilename);

void BuildContigAdjacency3(hashtable3 *ht, struct contigs_info *contigs_info,int K_size, string ContigFilename);

void BuildContigAdjacency4(hashtable4 *ht, struct contigs_info *contigs_info,int K_size, string ContigFilename);

void BuildContigAdjacency0(hashtable0 *ht, struct contigs_info *contigs_info,int K_size, string ContigFilename);

int BFSearchDist(struct hashtable* ht,struct hashtable* merge_ht, struct bucket* bktptr,struct bucket* obj_bktptr,int K_size, stacked_bucket &kmer_stack_beg,int max_depth,int max_dist);

int BFSearchDist2(struct hashtable2* ht,struct hashtable2* merge_ht, struct bucket2* bktptr,struct bucket2* obj_bktptr,int K_size, stacked_bucket2 &kmer_stack_beg,int max_depth,int max_dist);

int BFSearchDist3(struct hashtable3* ht,struct hashtable3* merge_ht, struct bucket3* bktptr,struct bucket3* obj_bktptr,int K_size, stacked_bucket3 &kmer_stack_beg,int max_depth,int max_dist);

int BFSearchDist4(struct hashtable4* ht,struct hashtable4* merge_ht, struct bucket4* bktptr,struct bucket4* obj_bktptr,int K_size, stacked_bucket4 &kmer_stack_beg,int max_depth,int max_dist);

void AppendMergeHT(hashtable *ht,hashtable *merge_ht);

void AppendMergeHT2(hashtable2 *ht,hashtable2 *merge_ht);

void AppendMergeHT3(hashtable3 *ht,hashtable3 *merge_ht);

void AppendMergeHT4(hashtable4 *ht,hashtable4 *merge_ht);

void AppendMergeHT0(hashtable0 *ht,hashtable0 *merge_ht,int Kmer_arr_sz);

void RemoveUnmappedNodes(hashtable *ht,hashtable2 *ht2,int K_size);

void RemoveUnmappedNodes3(hashtable3 *ht,int K_size);

void RemoveUnmappedNodes4(hashtable4 *ht,int K_size);

void ContigGapEst(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<int> &insert_sz_vt,vector<string>& filenames_vt,vector<bool> &LongLib,struct contigs_info * contigs_info,string ContigFilename,bool ResumePE,int64_t totReads,bool MatePair);

void ContigGapEst3(struct hashtable3 *ht,struct hashtable3 *merge_ht,int K_size,vector<int> &insert_sz_vt,vector<string>& filenames_vt,vector<bool> &LongLib,struct contigs_info * contigs_info,string ContigFilename,bool ResumePE,int64_t totReads,bool MatePair);

void ContigGapEst4(struct hashtable4 *ht,struct hashtable4 *merge_ht,int K_size,vector<int> &insert_sz_vt,vector<string>& filenames_vt,vector<bool> &LongLib,struct contigs_info * contigs_info,string ContigFilename,bool ResumePE,int64_t totReads,bool MatePair);

bool BackCheckLoop_ctg(int new_ctg,int end_node, map<int,struct BFS_path_info_ctg > & Visited_Path );

void BFSearchPathFinder(struct contigs_info *contigs_info,list<int> ctg_stack,map<int,list<int> > &dist_ctg,map<int,int > & unitig_dist,int dist_searched,map<int,vector<int> > &node_cov,int max_depth,bool GapClosingMode,bool MatePair);

void ResolvingRepeatsPE(vector<int> &insert_sz_vt,vector<string>& filenames_vt,struct contigs_info * contigs_info,string ContigFilename,int LinkCovTh,int UniqueLenTh,int ExpCov);

void ConstructContigGraph(struct hashtable *ht1,struct hashtable *merge_ht1, int K_size,contigs_info * contigs_info,string ContigFilename);

void ConstructContigGraph0(struct hashtable0 *ht,struct hashtable0 *merge_ht, int K_size, contigs_info * contigs_info,string ContigFilename);

#endif
