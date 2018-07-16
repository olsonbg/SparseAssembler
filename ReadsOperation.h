#ifndef __READS_OPERATION_H
#define __READS_OPERATION_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include "BasicDataStructure.h"

extern double UPDATERATE;

using namespace std;


struct SuperRead_t
{
	int PathsCnt;
	uint64_t idx;
	int PathsLength;
	int WellID;
	char seq[10000];
	bool PathFound;
	uint16_t depth;
	int ListSize;
	string extension;

};


struct SuperRead_ctg
{
	int PathsCnt;
	uint64_t idx;
	int PathsLength;

	bool PathFound;
	uint16_t depth;
	vector< vector<int> > Paths;

};


struct search_info
{
	int pos1,pos2;
	bool RightSearch;
	bool Flip_End;
	int offset;
};


struct BFS_path_info_v2
{
	int cov;
	int depth;
	int len;
	bool BothEndsUsed;
	uint64_t last_kmer;
	struct edge_node* last_bkt_edge;
};
struct BFS_path_info_v3
{
	int cov;
	int depth;
	int len;
	//bool BothEndsUsed;
	int RepeatVisits;
	int last_ctg;

};

struct reads_overlap_info
{
	vector< map<int32_t, vector<int32_t> > > left_overlaps,right_overlaps;
	vector<int> cov_vt;
	vector<bool> contained_vt,used_vt;

};

struct KmerInContig
{
	uint32_t contig_no : 31, flip : 1;
	int pos;
};
struct LongReadContigIndex
{
	map<int, KmerInContig> LR2CTG;
	map<int, vector<int> > CTG2LR, CTG2LR_2;
	map<int, vector<int> > layout;
	int KmerCovTh;
	int nMatches;
	//	vector<int> matching_positions_LR;
	//vector<KmerInContig> matching_positions_ctg;
};


struct BFS_reads_info
{
	int cov;
	int depth;
	int len;
	int last_read;
	vector<int> edge;
};

void SavingReadsTable(reads_table *reads_table);

void LoadingReadsTable(reads_table *reads_table);

void CollectingNonContainedReads(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename);


void CollectingNonContainedReads0(struct hashtable0 *ht0,struct hashtable0 *merge_ht0, int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename);


int BFSearchGapCloser_ctg(struct contigs_info *contigs_info,int max_depth,int max_dist,vector<int> &ctgs_in_current_read,SuperRead_ctg *SuperRead);

void ConstuctRefinedContigGraph(struct hashtable *ht1, int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename);

void ConstuctRefinedContigGraph0(struct hashtable0 *ht, int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename);


void CollectingNonContainedReadsSlow(struct hashtable *ht1, struct hashtable2 *ht2, int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename);

void CollectingNonContainedReadsSlow0(struct hashtable0 *ht, int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename);

void CollectingNonContainedPairsSlow(struct hashtable *ht1, struct hashtable2 *ht2, int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename);

void CollectingNonContainedPairsSlow0(struct hashtable0 *ht, int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename);


bool isSimplePath_read(reads_overlap_info *reads_overlap_info,int current_read,map<int,struct BFS_reads_info > &Visited_Path , map<int, int > &stacked_nodes);

void BreakLinks_read( reads_overlap_info *reads_overlap_info, map<int, int > &stacked_nodes,int node1, int node2);

void BacktrackBubbleRemoval_read(reads_overlap_info *reads_overlap_info,int last_read,int beg_read,map<int,struct BFS_reads_info > & Visited_Path , map<int ,int > &stacked_nodes);

void BFSearchBubbleRemoval_read(reads_overlap_info *reads_overlap_info,int beg_read,int max_depth);

void ConstructReadsOverlaps(string reads_info_name);

int BFSearchGapCloser_v2(struct hashtable* ht, uint64_t beg_kmer,uint64_t end_kmer,int K_size,int max_depth,int max_dist,SuperRead_t *SuperRead);

void CollectingNonContainedPairs(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<string>& filenames_vt,struct contigs_info * contigs_info,string ContigFilename);

#endif
