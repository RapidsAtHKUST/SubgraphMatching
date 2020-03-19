#ifndef _INTER_ALGOS_H
#define _INTER_ALGOS_H

#include "utils/util.hpp"

// ScalarMerge:
int intersect_scalarmerge_uint(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
// ScalarMerge+BSRSet:
int intersect_scalarmerge_bsr(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c);

// ScalarGalloping:
int intersect_scalargalloping_uint(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
// ScalarGalloping+BSRSet:
int intersect_scalargalloping_bsr(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c);

// SIMDGalloping:
int intersect_simdgalloping_uint(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
// SIMDGalloping+BSRSet:
int intersect_simdgalloping_bsr(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c);

// QFilter:
int intersect_qfilter_uint_b4(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
int intersect_qfilter_uint_b4_v2(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c);

// QFilter+BSRSet:
int intersect_qfilter_bsr_b4(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c);
int intersect_qfilter_bsr_b4_v2(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c);

// Shuffling:
int intersect_shuffle_uint_b4(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
int intersect_shuffle_uint_b8(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);

// Shuffling+BSRSet:
int intersect_shuffle_bsr_b4(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c);

// BMiss:
int intersect_bmiss_uint_b4(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
// BMiss+STTNI (block size = 8):
int intersect_bmiss_uint_sttni_b8(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
// HieraInter:
int intersect_hierainter_uint_sttni(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c);
int hierainter_offline_partition(int *set_a, int size_a, uint16_t *hi_set);
int hierainter_offline_combine(uint16_t *hi_set, int hi_size, int *set_a);
int hierainter_online_intersect_low16bit(uint16_t *set_a, int size_a,
            uint16_t *set_b, int size_b, uint16_t *set_c);
int hierainter_online_intersect_high16bit(uint16_t *par_a, int size_a,
            uint16_t *par_b, int size_b, uint16_t *par_c);

int offline_uint_trans_bsr(int *set_a, int size_a, int *bases_a, int *states_a);
int offline_bsr_trans_uint(int *bases_a, int *states_a, int size_a, int *set_a);
#endif

