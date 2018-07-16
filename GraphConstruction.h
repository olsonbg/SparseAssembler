#ifndef __GRAPH_CONSTRUCTION_H
#define __GRAPH_CONSTRUCTION_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include "BasicDataStructure.h"

// initialize a hashtable
void Init_HT(struct hashtable* ht,size_t ht_sz);

void Init_HT2(struct hashtable2* ht,size_t ht_sz);

void Init_HT3(struct hashtable3* ht,size_t ht_sz);

void Init_HT4(struct hashtable4* ht,size_t ht_sz);

void Init_HT0(struct hashtable0* ht,size_t ht_sz);

//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 2
bool look_up_in_a_list(uint64_t seq,struct bucket *** ptr);

bool look_up_in_a_list2(struct kmer_t2 *seq,struct bucket2 *** ptr);

bool look_up_in_a_list3(struct kmer_t3 *seq,struct bucket3 *** ptr);

bool look_up_in_a_list4(struct kmer_t4 *seq,struct bucket4 *** ptr);


bool look_up_in_a_list0( uint64_t *seq,struct bucket0 *** ptr,int arr_sz);

bool look_up_in_a_list_rm(uint64_t seq,struct bucket_rm *** ptr);

bool look_up_in_a_list_rm2(struct kmer_t2 *seq,struct bucket_rm2 *** ptr);

bool look_up_in_a_list_rm3(struct kmer_t3 *seq,struct bucket_rm3 *** ptr);

bool look_up_in_a_list_rm4(struct kmer_t4 *seq,struct bucket_rm4 *** ptr);

bool look_up_in_a_list_rm0(uint64_t *seq,struct bucket_rm0 *** ptr,int Kmer_arr_sz);

// compare 2 arrays.
int uint64_t_cmp(uint64_t* A,uint64_t* B,int Kmer_arr_sz);

//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 1
bool look_up_in_a_list_r1(uint64_t seq,struct bucket_r1 *** ptr);

bool look_up_in_a_list2_r1(struct kmer_t2 *seq,struct bucket2_r1 *** ptr);

bool look_up_in_a_list3_r1(struct kmer_t3 *seq,struct bucket3_r1 *** ptr);

bool look_up_in_a_list4_r1(struct kmer_t4 *seq,struct bucket4_r1 *** ptr);

bool look_up_in_a_list0_r1(uint64_t *seq,struct bucket0_r1 *** ptr,int arr_sz);

// sparse k-mer graph construction function, it contains 2 rounds. In the first round we select the sparse k-mers,
//in the second, we build the links between the k-mers.
void Sparse_Kmer_Graph_Construction(struct read_t *read,struct hashtable *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info ,int round);

void Sparse_Kmer_Graph_Construction0(struct read_t *read,struct hashtable0 *ht, key_table *key_table, int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info ,int round);

void Sparse_Kmer_Graph_Construction2(struct read_t *read,struct hashtable2 *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info,int round);

void Sparse_Kmer_Graph_Construction3(struct read_t *read,struct hashtable3 *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info, int round);

void Sparse_Kmer_Graph_Construction4(struct read_t *read,struct hashtable4 *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info,int round);

void Sparse_Kmer_Index_Construction(struct read_t *read,struct hashtable *ht,int64_t *bucket_count,int K_size,int gap,BF_info * BF_info ,int round,struct read_index *read_index);

//convert the bucket type from round 1 to round 2. The buckets in round 2 are more expensive
void SwitchBuckets(hashtable *ht,hashtable2 *ht2,int K_size);
void SwitchBuckets3(hashtable3 *ht3,int K_size);


void SwitchBuckets4(hashtable4 *ht4,int K_size);

void SwitchBuckets0(hashtable0 *ht0);

void OutputSparseKmers(hashtable *ht,int K_size,bool Bloom);

void OutputSparseKmers2(hashtable2 *ht,int K_size,bool Bloom);

void OutputSparseKmers3(hashtable3 *ht,int K_size,bool Bloom);
void OutputSparseKmers4(hashtable4 *ht,int K_size,bool Bloom);

//save the sparse graph into disk
void SavingSparseKmerGraph(hashtable *ht,string &fname);

//save the bubble removal information into disk
void SavingMergeHT(hashtable *ht);

void SavingSparseKmerGraph_E2(hashtable *ht,string &fname);

void SavingSparseKmerGraph2(hashtable2 *ht,string &fname);

void SavingMergeHT2(hashtable2 *ht);

void SavingSparseKmerGraph3(hashtable3 *ht,string &fname);

void SavingMergeHT3(hashtable3 *ht);

void SavingSparseKmerGraph4(hashtable4 *ht,string &fname);

void SavingMergeHT4(hashtable4 *ht);

void SavingSparseKmerGraph0(hashtable0 *ht,string &fname,int Kmer_arr_sz);

void SavingMergeHT0(hashtable0 *ht,int Kmer_arr_sz);

void SavingSparseKmerIndex(hashtable *ht,read_index * read_index, string &fname);

void LoadingSparseKmerIndex(hashtable *ht,read_index * read_index, string &fname);

//load the saved information
void LoadingSparseKmerGraph(hashtable *ht,string &fname);

void LoadingMergeHT(hashtable *ht);

void LoadingSparseKmerGraph2(hashtable2 *ht,string &fname);

void LoadingMergeHT2(hashtable2 *ht);

void LoadingSparseKmerGraph3(hashtable3 *ht,string &fname);

void LoadingMergeHT3(hashtable3 *ht);

void LoadingSparseKmerGraph4(hashtable4 *ht,string &fname);

void LoadingMergeHT4(hashtable4 *ht);


void LoadingSparseKmerGraph0(hashtable0 *ht,key_table *key_table,string &fname,int Kmer_arr_sz);

void LoadingMergeHT0(hashtable0 *ht,key_table *key_table,int Kmer_arr_sz);

// free the graph
void FreeSparseKmerGraph(struct hashtable *ht);

void FreeSparseKmerGraph2(struct hashtable2 *ht);

void FreeSparseKmerGraph3(struct hashtable3 *ht);

void FreeSparseKmerGraph4(struct hashtable4 *ht);

void FreeSparseKmerGraph0(struct hashtable0 *ht,key_table *key_table);


void ScanDataset(vector<string > & in_filenames_vt,uint64_t *tot_bases,uint64_t *numReads,uint64_t *totReads,int *maxlen);

#endif
