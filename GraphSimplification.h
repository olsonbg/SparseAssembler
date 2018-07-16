#ifndef __GRAPH_SIMPLIFICATION_H
#define __GRAPH_SIMPLIFICATION_H
#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include "BasicDataStructure.h"
#include "GraphConstruction.h"

//break the bubble links in the bfs.

void BreakLinks( map<struct bucket* ,int > &stacked_nodes, struct bucket* bktptr1, struct bucket* bktptr2,int K_size,int edge_len);

void BreakLinks2( map<struct bucket2* ,int > &stacked_nodes, struct bucket2* bktptr1, struct bucket2* bktptr2,int K_size,int edge_len);

void BreakLinks3( map<struct bucket3* ,int > &stacked_nodes, struct bucket3* bktptr1, struct bucket3* bktptr2,int K_size,int edge_len);


void BreakLinks4( map<struct bucket4* ,int > &stacked_nodes, struct bucket4* bktptr1, struct bucket4* bktptr2,int K_size,int edge_len);

void BreakLinks0( map<struct bucket0* ,int > &stacked_nodes, struct bucket0* bktptr1, struct bucket0* bktptr2,int K_size,int edge_len);

void RemovingWeakNodesAndEdges(hashtable *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt);

void RemovingWeakNodesAndEdges2(hashtable2 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt);


void RemovingWeakNodesAndEdges3(hashtable3 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt);


void RemovingWeakNodesAndEdges4(hashtable4 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt);

void RemovingWeakNodesAndEdges0(hashtable0 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt);

void RemovingWeakNodes_r1(hashtable *ht,hashtable2 *ht2,int K_size,int NodeCovTh, int64_t *bucket_cnt);

void RemovingWeakNodes3_r1(hashtable3 *ht3,int K_size,int NodeCovTh, int64_t *bucket_cnt);

void RemovingWeakNodes4_r1(hashtable4 *ht4,int K_size,int NodeCovTh, int64_t *bucket_cnt);

void RemovingWeakNodes0_r1(hashtable0 *ht0, int NodeCovTh, int64_t *bucket_cnt);

void MergeNode(hashtable *merge_ht,uint64_t kmer,uint64_t merge_kmer,bool flip);

void MergeNode2(hashtable2 *merge_ht,kmer_t2 kmer,kmer_t2 merge_kmer,bool flip);

void MergeNode3(hashtable3 *merge_ht,kmer_t3 kmer,kmer_t3 merge_kmer,bool flip);

void MergeNode4(hashtable4 *merge_ht,kmer_t4 kmer,kmer_t4 merge_kmer,bool flip);

void MergeNode0(hashtable0 *merge_ht,uint64_t *kmer,uint64_t *merge_kmer,bool flip,int Kmer_arr_sz);

bool isSimplePath(struct bucket* bktptr,map<bucket*,struct BFS_path_info > & Visited_Path , map<struct bucket* ,int > &stacked_nodes);

bool isSimplePath2(struct bucket2* bktptr,map<bucket2*,struct BFS_path_info2 > & Visited_Path , map<struct bucket2* ,int > &stacked_nodes);

bool isSimplePath3(struct bucket3* bktptr,map<bucket3*,struct BFS_path_info3 > & Visited_Path , map<struct bucket3* ,int > &stacked_nodes);

bool isSimplePath4(struct bucket4* bktptr,map<bucket4*,struct BFS_path_info4 > & Visited_Path , map<struct bucket4* ,int > &stacked_nodes);

bool isSimplePath0(hashtable0 *ht,struct bucket0* bktptr,map<bucket0*,struct BFS_path_info0 > & Visited_Path , map<struct bucket0* ,int > &stacked_nodes);

// backtrack operation in the bubble removal
void BacktrackBubbleRemoval(struct hashtable* merge_ht,struct bucket* bktptr,struct bucket* bktptr_merged,struct bucket* beg_bkt,map<bucket*,struct BFS_path_info > & Visited_Path , map<struct bucket* ,int > &stacked_nodes, int K_size);

//check for self loops in the bubble removal
bool BackCheckLoop(struct bucket* bktptr,struct bucket* end_bkt,map<bucket*,struct BFS_path_info > & Visited_Path );

void BFSearchBubbleRemoval(struct hashtable* ht,struct hashtable* merge_ht, struct bucket* bktptr,int K_size, list<struct stacked_bucket>& kmer_stack,int PathCovTh,int max_depth,int PathSim);

void BacktrackBubbleRemoval2(struct hashtable2* merge_ht,struct bucket2* bktptr,struct bucket2* bktptr_merged,struct bucket2* beg_bkt,map<bucket2*,struct BFS_path_info2 > & Visited_Path , map<struct bucket2* ,int > &stacked_nodes, int K_size);

bool BackCheckLoop2(struct bucket2* bktptr,struct bucket2* end_bkt,map<bucket2*,struct BFS_path_info2 > & Visited_Path );

// search operation in the bubble removal
void BFSearchBubbleRemoval2(struct hashtable2* ht,struct hashtable2* merge_ht,struct bucket2* bktptr,int K_size, list<struct stacked_bucket2>& kmer_stack,int PathCovTh,int max_depth,int PathSim);

void BacktrackBubbleRemoval3(struct hashtable3* merge_ht, struct bucket3* bktptr,struct bucket3* bktptr_merged,struct bucket3* beg_bkt,map<bucket3*,struct BFS_path_info3 > & Visited_Path , map<struct bucket3* ,int > &stacked_nodes, int K_size);

bool BackCheckLoop3(struct bucket3* bktptr,struct bucket3* end_bkt,map<bucket3*,struct BFS_path_info3 > & Visited_Path );

void BFSearchBubbleRemoval3(struct hashtable3* ht,struct hashtable3* merge_ht,struct bucket3* bktptr,int K_size, list<struct stacked_bucket3>& kmer_stack,int PathCovTh,int max_depth,int PathSim);

void BacktrackBubbleRemoval4(struct hashtable4* merge_ht,struct bucket4* bktptr,struct bucket4* bktptr_merged,struct bucket4* beg_bkt,map<bucket4*,struct BFS_path_info4 > & Visited_Path , map<struct bucket4* ,int > &stacked_nodes, int K_size);

bool BackCheckLoop4(struct bucket4* bktptr,struct bucket4* end_bkt,map<bucket4*,struct BFS_path_info4 > & Visited_Path );

void BFSearchBubbleRemoval4(struct hashtable4* ht,struct hashtable4* merge_ht,struct bucket4* bktptr,int K_size,list<struct stacked_bucket4>& kmer_stack,int PathCovTh,int max_depth,int PathSim);

void BacktrackBubbleRemoval0(struct hashtable0* merge_ht,struct bucket0* bktptr,struct bucket0* bktptr_merged,struct bucket0* beg_bkt,map<bucket0*,struct BFS_path_info0 > & Visited_Path , map<struct bucket0* ,int > &stacked_nodes, int K_size);

bool BackCheckLoop0(struct bucket0* bktptr,struct bucket0* end_bkt,map<bucket0*,struct BFS_path_info0 > & Visited_Path );

void BFSearchBubbleRemoval0(struct hashtable0* ht,struct hashtable0* merge_ht,struct bucket0* bktptr,int K_size, list<struct stacked_bucket0>& kmer_stack,int PathCovTh,int max_depth,int PathSim);

// call the bubble removal
void GraphSimplification(hashtable *ht,hashtable * merge_ht, hashtable2 *ht2, hashtable2 * merge_ht2,int K_size,int PathCovTh,int max_depth,int PathSim);

void GraphSimplification3(hashtable3 *ht3, hashtable3 * merge_ht3,int K_size,int PathCovTh,int max_depth,int PathSim);

void GraphSimplification4(hashtable4 *ht4, hashtable4 * merge_ht4,int K_size,int PathCovTh,int max_depth,int PathSim);

void GraphSimplification0(hashtable0 *ht, hashtable0 * merge_ht,int K_size,int PathCovTh,int max_depth,int PathSim);

#endif
