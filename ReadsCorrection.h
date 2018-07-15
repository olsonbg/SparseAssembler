#ifndef __ReadsCorrection_H
#define __ReadsCorrection_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include "BasicDataStructure.h"

inline void replace_one_nc(uint64_t * bitsarr_in,int bitsarr_len,int begin_pos,uint64_t nt);

int BFsearch_denoising(struct read_t *read,struct hashtable *ht,struct hashtable2 *ht2,size_t *CovTh,size_t *CorrTh,int K_size,int gap,int *last_correct,bool SEARCH_RIGHT);

bool Sparse_Denoising(struct read_t *read,struct hashtable *ht,struct hashtable2 *ht2,size_t *CovTh,size_t *CorrTh,int K_size,int gap,uint64_t *correction_cnt,bool Hybrid);

#endif
