#ifndef __BASIC_DATA_STRUCTURE_H
#define __BASIC_DATA_STRUCTURE_H

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include "OpenFile.h"

using namespace std;
//typedef unsigned __int64 uint64_t;
//typedef __int64 int64_t;
//typedef unsigned __int8 uint8_t;
//typedef unsigned __int16 uint16_t;
//typedef unsigned __int32 uint32_t;
//typedef __int32 int32_t;



// These are the structures to save the k-mers.32 bases, 64 bases, 96 bases, 128 bases.
struct kmer_t
{
	uint64_t kmer;
};

struct kmer_t2
{
	uint64_t kmer[2];
};

struct kmer_t3
{
	uint64_t kmer[3];
};

struct kmer_t4
{
	uint64_t kmer[4];
};




// edge structure, maximum length 25
struct edge_node
{
	uint64_t edge:50,edge_cov:4,len:6,masked:1,removed:1;
	struct edge_node *nxt_edge;
};
// edge structure, maximum length 64
struct edge_node2
{
	uint64_t edge[2];
	uint16_t edge_cov:8,len:8;
	struct edge_node2 *nxt_edge;
};

// information for a k-mer node
struct kmer_info
{
	//flags, recording if the node is used in a search, whether the node have branches on its sides...
	uint8_t used:1,split_left:1,split_right:1,removed:1,flip:1,marked:1,repeat:1,masked:1;
	//pointers to edge links
	struct edge_node *left;
	struct edge_node *right;
	uint16_t cov1;
	//for building scaffolds
	int32_t contig_no;
	int32_t cod;
};

//split the above structure, round 1 of graph construction:
struct kmer_info_r1
{
	uint16_t cov1;
};


//Bloom filter structure
struct BF_info
{
	uint8_t * BF_HT;
	uint64_t m;
	int d;
	bool Bloom;
};

// a complete k-mer node,
struct bucket
{
	struct kmer_t kmer_t;	//32 bp
	struct kmer_info kmer_info;
	bucket *nxt_bucket;
};

struct bucket2
{
	struct kmer_t2 kmer_t2;	//64 bp
	struct kmer_info kmer_info;
	bucket2 *nxt_bucket;
};

struct bucket3
{
	struct kmer_t3 kmer_t3;
	struct kmer_info kmer_info;
	bucket3 *nxt_bucket;
};

struct bucket4
{
	struct kmer_t4 kmer_t4;
	struct kmer_info kmer_info;
	bucket4 *nxt_bucket;
};


struct bucket0
{
	uint64_t *kmer_t;	//any kmer size
	struct kmer_info kmer_info;
	bucket0 *nxt_bucket;
};




// a structure recording the information of removed k-mers in the BFS bubble removal

struct bucket_rm
{
	struct kmer_t kmer_t;
	struct kmer_t merged_kmer;
	bool flip;
	bucket_rm *nxt_bucket;
};

struct bucket_rm2
{
	struct kmer_t2 kmer_t2;
	struct kmer_t2 merged_kmer;
	bool flip;
	bucket_rm2 *nxt_bucket;
};
struct bucket_rm3
{
	struct kmer_t3 kmer_t3;
	struct kmer_t3 merged_kmer;
	bool flip;
	bucket_rm3 *nxt_bucket;
};
struct bucket_rm4
{
	struct kmer_t4 kmer_t4;
	struct kmer_t4 merged_kmer;
	bool flip;
	bucket_rm4 *nxt_bucket;
};

struct bucket_rm0
{
	uint64_t *kmer_t;
	uint64_t *merged_kmer;
	bool flip;
	bucket_rm0 *nxt_bucket;
};

//bucket in round 1
struct bucket_r1
{
	struct kmer_t kmer_t;
	struct kmer_info_r1 kmer_info;
	bucket_r1 *nxt_bucket;
};

struct bucket2_r1
{
	struct kmer_t2 kmer_t2;
	struct kmer_info_r1 kmer_info;
	bucket2_r1 *nxt_bucket;
};

struct bucket3_r1
{
	struct kmer_t3 kmer_t3;
	struct kmer_info_r1 kmer_info;
	bucket3_r1 *nxt_bucket;
};

struct bucket4_r1
{
	struct kmer_t4 kmer_t4;
	struct kmer_info_r1 kmer_info;
	bucket4_r1 *nxt_bucket;
};

struct bucket0_r1
{
	uint64_t * kmer_t;
	struct kmer_info_r1 kmer_info;
	bucket0_r1 *nxt_bucket;
};
//hashtable structure

struct hashtable
{
	struct bucket **store_pos;
	size_t ht_sz;
};

struct hashtable2
{
	struct bucket2 **store_pos;
	size_t ht_sz;
};


struct hashtable3
{
	struct bucket3 **store_pos;
	size_t ht_sz;
};

struct hashtable4
{
	struct bucket4 **store_pos;
	size_t ht_sz;
};

struct hashtable0
{
	struct bucket0 **store_pos;
	size_t ht_sz;
};

//read structure

struct read_t
{
	char tag[1000];
	bool error_nt[1000];
	char c_seq[10000];//char representation
	//char *c_seq;//char representation

	//uint64_t read_bits[10000];//bit representation
	uint64_t *read_bits;
	//char read[1000];//char representation
	int readLen;// read length
	int read_idx;
};

struct key_table
{
	bool in_use;
	list<uint64_t *> pblocks;
	int BytesPerKey;
	int KeysPerBlock;
	int current_block;
	int current_index;
};






struct read_index
{
	vector< map<uint64_t,bool> > repeat_maps;
	uint64_t repeat_cnt;
	int MaxMatch;
	vector<int> read_len_vt;
};


struct reads_table
{
	bool in_use;
	list<uint64_t *> pblocks;
	int BytesPerBlock;
	int current_block;
	int current_byte;
	int current_read;
	map<int,uint64_t *> read_ptr;
	vector<int> read_len_vt;
};

//contig graph

struct contigs_info
{
	int total_contigs;
	int K_size;
	vector<int> contig_sz_vt,kmer_cnt_vt,comp_vt;
	vector<int> contigs_hp_b,contigs_hp_e;
	vector<string> contigs_str;
	map<int, vector<int> > Length_ID;
	map<int, vector<int> > Cov_Length;
	//map<int, vector<int> > ctg_in_scf;
	//vector<vector<int> > scaffolds;
	//vector<vector<int> > gaps_in_scaffolds;

	vector < vector<int>::iterator> LengthRank;
	vector<int> cov_vt;
	vector<struct c_info> c_info_vt;
	vector< map<int,struct scaffold_contig_info> > scaffold_adjacency_left,scaffold_adjacency_right;
	vector< map<int,struct adjacent_contig_info> > contig_adjacency_left,contig_adjacency_right;
};




//path information in the BFS bubble removal

struct BFS_path_info
{
	int cov;
	int depth;
	int len;
	bool BothEndsUsed;
	struct bucket* last_bkt;
	struct edge_node* last_bkt_edge;
};





struct BFS_path_info2
{
	int cov;
	int depth;
	int len;
	struct bucket2* last_bkt;
	struct edge_node* last_bkt_edge;
};

struct BFS_path_info3
{
	int cov;
	int depth;
	int len;
	struct bucket3* last_bkt;
	struct edge_node* last_bkt_edge;
};

struct BFS_path_info4
{
	int cov;
	int depth;
	int len;
	struct bucket4* last_bkt;
	struct edge_node* last_bkt_edge;
};

struct BFS_path_info0
{
	int cov;
	int depth;
	int len;
	struct bucket0* last_bkt;
	struct edge_node* last_bkt_edge;
};


//outdated

struct BFS_path_info_ctg
{
	int cov;
	int depth;
	int len;
	int last_ctg;
};


//a stack in the BFS
struct stacked_bucket
{
	bucket *bktptr;
	bool RightSearch;
	bool BothEndsUsed;
};
struct stacked_bucket2
{
	bucket2 *bktptr;
	bool RightSearch;
};
struct stacked_bucket3
{
	bucket3 *bktptr;
	bool RightSearch;
};
struct stacked_bucket4
{
	bucket4 *bktptr;
	bool RightSearch;
};

struct stacked_bucket0
{
	bucket0 *bktptr;
	bool RightSearch;
};


bool get_a_fasta_read(ifstream & fasta_in, string &tag, string &str, string & n_tag);

bool get_a_fasta_read(boost::iostreams::filtering_stream<boost::iostreams::input> *fasta_in, string &tag, string &str, string & n_tag);


bool get_a_fastq_read(ifstream & fastq_in, string &tag, string &seq, string & quality);

bool get_a_fastq_read(boost::iostreams::filtering_stream<boost::iostreams::input> *fastq_in, string & tag, string &seq, string &quality);


bool basic_quality_check(string &seq_s);




//left shift and right shift of shift_sz bits of the whole bit array, arr_sz is the array length
static inline void L_shift_NB(uint64_t * bitsarr, int shift_sz,int arr_sz)
{

	uint64_t temp_arr[100];

	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/

	int jmp=shift_sz/64;
	int offset=shift_sz%64;

	for (int i=0;i<arr_sz;++i)
	{
		if(i+jmp+1<arr_sz)
		{

			uint64_t tt=0;
			if(offset==0)
			{
				tt=0;
			}
			else
			{
				tt=(bitsarr[i+jmp+1]>>(64-offset));
			}
			temp_arr[i]=((bitsarr[i+jmp]<<offset)|tt);
		}
		if(i+jmp+1==arr_sz)
		{temp_arr[i]=bitsarr[i+jmp]<<offset;}
		if(i+jmp+1>arr_sz)
		{temp_arr[i]=0;}

	}

	memcpy(bitsarr,temp_arr,sizeof(uint64_t)*arr_sz);
/*
	for (int i=0;i<arr_sz;++i)
	{
		bitsarr[i]=temp_arr[i];
	}
	*/
}

static inline void R_shift_NB(uint64_t * bitsarr, int shift_sz,int arr_sz)
{
	uint64_t temp_arr[100];
	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/
	int jmp=shift_sz/64;
	int offset=shift_sz%64;

	for (int i=arr_sz-1;i>=0;--i)
	{
		if(i-jmp>0)
		{
			if(offset>0)
			{temp_arr[i]=(bitsarr[i-jmp]>>offset)|(bitsarr[i-jmp-1]<<(64-offset));}
			else
			{temp_arr[i]=bitsarr[i-jmp];}
		}
		if (i-jmp==0)
		{
			if(offset>0)
			{temp_arr[i]=(bitsarr[i-jmp]>>offset);}
			else
			{temp_arr[i]=bitsarr[i-jmp];}
		}
		if (i-jmp<0)
		{temp_arr[i]=0;}

	}
	memcpy(bitsarr,temp_arr,sizeof(uint64_t)*arr_sz);
	/*
	for (int i=0;i<arr_sz;++i)
	{
		bitsarr[i]=temp_arr[i];
	}
	*/

}

// get reverse complement of a k-mer.
static inline uint64_t get_rev_comp_seq(uint64_t seq, int seq_size)
{
	seq =~seq;

	seq = ((seq & 0x3333333333333333 )<< 2) | ((seq & 0xCCCCCCCCCCCCCCCC )>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0F )<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0 )>> 4);
	seq = ((seq & 0x00FF00FF00FF00FF )<< 8) | ((seq & 0xFF00FF00FF00FF00 )>> 8);
	seq = ((seq & 0x0000FFFF0000FFFF )<<16) | ((seq & 0xFFFF0000FFFF0000 )>>16);
	seq = ((seq & 0x00000000FFFFFFFF )<<32) | ((seq & 0xFFFFFFFF00000000 )>>32);

	return seq >> (64 - (seq_size*2));
}

static inline uint64_t* get_rev_comp_seq_arr(uint64_t *seq_arr, int seq_size,int arr_sz)
{


	int tot_bits=arr_sz*64;
	for(int i=0;i<arr_sz;++i)
	{
		seq_arr[i]=~seq_arr[i];
		seq_arr[i] = ((seq_arr[i] & 0x3333333333333333 )<< 2) | ((seq_arr[i] & 0xCCCCCCCCCCCCCCCC )>> 2);
		seq_arr[i] = ((seq_arr[i] & 0x0F0F0F0F0F0F0F0F )<< 4) | ((seq_arr[i] & 0xF0F0F0F0F0F0F0F0 )>> 4);
		seq_arr[i] = ((seq_arr[i] & 0x00FF00FF00FF00FF )<< 8) | ((seq_arr[i] & 0xFF00FF00FF00FF00 )>> 8);
		seq_arr[i] = ((seq_arr[i] & 0x0000FFFF0000FFFF )<<16) | ((seq_arr[i] & 0xFFFF0000FFFF0000 )>>16);
		seq_arr[i] = ((seq_arr[i] & 0x00000000FFFFFFFF )<<32) | ((seq_arr[i] & 0xFFFFFFFF00000000 )>>32);
	}

	int j=0,k=arr_sz-1;
	for (;j<k;++j,--k)
	{
		uint64_t temp;
		temp=seq_arr[j];
		seq_arr[j]=seq_arr[k];
		seq_arr[k]=temp;
	}

	R_shift_NB(seq_arr,tot_bits-(seq_size*2),arr_sz);
	return seq_arr;
	//return seq >> (64 - (seq_size*2));
}

// get sub bit array from a bit array.
 inline void get_sub_arr(uint64_t * bitsarr_in,int bitsarr_len,int begin_pos,int sub_sz,uint64_t * bitsarr_out)
{
	if(bitsarr_len<sub_sz)
	{cout<<"Error! Input kmer too short."<<bitsarr_len <<" "<<sub_sz<<endl;return;}
	int arr_sz_in=bitsarr_len/32+1;
	int rem=bitsarr_len%32;
	if(rem==0)
	{arr_sz_in--;}

	int arr_sz_out=sub_sz/32+1;
	if(sub_sz%32==0)
	{arr_sz_out--;}

	uint64_t temp_arr[10];
	memset(temp_arr,0,sizeof(temp_arr));

	memset(bitsarr_out,0,sizeof(uint64_t)*arr_sz_out);

	int rem2=(32-rem+begin_pos)%32;
	int block_beg=(32-rem+begin_pos)/32;
	if(rem==0)
	{block_beg--;}

	int rem3=(32-rem+begin_pos+sub_sz)%32;
	int block_end=(32-rem+begin_pos+sub_sz)/32;
	if(rem3!=0)
	{rem3=32-rem3;}
	else
	{
		block_end--;
	}
	if(rem==0)
	{block_end--;}

	int orig_sz=(block_end-block_beg+1);
	memcpy(temp_arr,&bitsarr_in[block_beg],orig_sz*sizeof(uint64_t));
	L_shift_NB(temp_arr,rem2*2,orig_sz);
	R_shift_NB(temp_arr,(rem2+rem3)%32*2,arr_sz_out);
	memcpy(bitsarr_out,temp_arr,arr_sz_out*sizeof(uint64_t));


}

uint64_t* str2bitsarr(const char * c_str,int len, uint64_t* b_str,int arr_sz );

 //hash functions
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed );

uint64_t MurmurHash64B ( const void * key, int len, unsigned int seed );

 //convert a string of nucleotide bases into bit array.
void Init_Read(string &seq,struct read_t & read);



//get the complement of a string of nucleotide bases
void complement_str(string & str);

//convert a bit array into a string
char * bitsarr2str(uint64_t* b_seq, int len,char * c_str,int arr_sz);

void string_shrinkage(string &s);

#endif
