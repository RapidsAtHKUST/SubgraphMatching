#include "intersection_algos.hpp"

#ifdef USE_ALIGNED
#define MM_LOAD_SI_128 _mm_load_si128
#define MM_STORE_SI_128 _mm_store_si128
#else
#define MM_LOAD_SI_128 _mm_loadu_si128
#define MM_STORE_SI_128 _mm_storeu_si128
#endif

constexpr int cyclic_shift1 = _MM_SHUFFLE(0,3,2,1); //rotating right
constexpr int cyclic_shift2 = _MM_SHUFFLE(2,1,0,3); //rotating left
constexpr int cyclic_shift3 = _MM_SHUFFLE(1,0,3,2); //between

static const __m128i all_zero_si128 = _mm_setzero_si128();
static const __m128i all_one_si128 = _mm_set_epi32(0xffffffff, 0xffffffff,
        0xffffffff, 0xffffffff);

static const uint8_t shuffle_pi8_array[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 
    8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 
    12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 
    8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
};
static const __m128i *shuffle_mask = (__m128i*)(shuffle_pi8_array);

uint32_t* prepare_shuffling_dict_avx()
{
    uint32_t* arr = new uint32_t[2048];
    for(int i=0; i<256; ++i){
        int count=0, rest=7;
        for(int b=0; b<8; ++b){
            if(i & (1 << b)){
                // n index at pos p - move nth element to pos p
                arr[i*8 + count] = b; // move all set bits to beginning
                ++count;
            }else{
                arr[i*8 + rest] = b; // move rest at the end
                --rest;
            }
        }
    }
    return arr;
}
static const uint32_t *shuffle_mask_avx =prepare_shuffling_dict_avx();

int intersect_scalarmerge_uint(int *set_a, int size_a, int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;    
}

int intersect_scalarmerge_bsr(int* bases_a, int* states_a, int size_a,
        int* bases_b, int* states_b, int size_b,
        int* bases_c, int* states_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_scalargalloping_uint(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    if (size_a > size_b)
        return intersect_scalargalloping_uint(set_b, size_b, set_a, size_a, set_c);
    
    int j = 0, size_c = 0;    
    for (int i = 0; i < size_a; ++i) {
        // double-jump:
        int r = 1;
        while (j + r < size_b && set_a[i] > set_b[j + r]) r <<= 1;
        // binary search:
        int right = (j + r < size_b) ? (j + r) : (size_b - 1);
        if (set_b[right] < set_a[i]) break;        
        int left = j + (r >> 1);
        while (left < right) {
            int mid = (left + right) >> 1;
            if (set_b[mid] >= set_a[i]) right = mid;
            else left = mid + 1;
        }
        j = left;

        if (set_a[i] == set_b[j]) set_c[size_c++] = set_a[i];
    }

    return size_c;
}

int intersect_scalargalloping_bsr(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c)
{
    if (size_a > size_b)
        return intersect_scalargalloping_bsr(bases_b, states_b, size_b,
                bases_a, states_a, size_a,
                bases_c, states_c);

    int j = 0, size_c = 0;
    for (int i = 0; i < size_a; ++i) {
        // double-jump:
        int r = 1;
        while (j + r < size_b && bases_a[i] > bases_b[j + r]) r <<= 1;
        // binary search:
        int right = (j + r < size_b) ? (j + r) : (size_b - 1);
        if (bases_b[right] < bases_a[i]) break;        
        int left = j + (r >> 1);
        while (left < right) {
            int mid = (left + right) >> 1;
            if (bases_b[mid] >= bases_a[i]) right = mid;
            else left = mid + 1;
        }
        j = left;

        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
        }
    }

    return size_c;
}

int intersect_simdgalloping_uint(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    if (size_a > size_b)
        return intersect_simdgalloping_uint(set_b, size_b, set_a, size_a, set_c);

    int i = 0, j = 0, size_c = 0;
    int qs_b = size_b - (size_b & 3);
    for (i = 0; i < size_a; ++i) {
        // double-jump:
        int r = 1;
        while (j + (r << 2) < qs_b && set_a[i] > set_b[j + (r << 2) + 3]) r <<= 1;
        // binary search:
        int upper = (j + (r << 2) < qs_b) ? (r) : ((qs_b - j - 4) >> 2);
        if (set_b[j + (upper << 2) + 3] < set_a[i]) break;        
        int lower = (r >> 1);
        while (lower < upper) {
            int mid = (lower + upper) >> 1;
            if (set_b[j + (mid << 2) + 3] >= set_a[i]) upper = mid;
            else lower = mid + 1;
        }
        j += (lower << 2);

        __m128i v_a = _mm_set_epi32(set_a[i], set_a[i], set_a[i], set_a[i]);
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));
        __m128i cmp_mask = _mm_cmpeq_epi32(v_a, v_b);
        int mask = _mm_movemask_ps((__m128)cmp_mask);
        if (mask != 0) set_c[size_c++] = set_a[i];
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_simdgalloping_bsr(int* bases_a, int* states_a, int size_a,
        int* bases_b, int* states_b, int size_b,
        int* bases_c, int* states_c)
{
    if (size_a > size_b)
        return intersect_simdgalloping_bsr(bases_b, states_b, size_b,
                bases_a, states_a, size_a,
                bases_c, states_c);

    int i = 0, j = 0, size_c = 0;
    int qs_b = size_b - (size_b & 3);
    for (i = 0; i < size_a; ++i) {
        // double-jump:
        int r = 1;
        while (j + (r << 2) < qs_b && bases_a[i] > bases_b[j + (r << 2) + 3]) r <<= 1;
        // binary search:
        int upper = (j + (r << 2) < qs_b) ? (r) : ((qs_b - j - 4) >> 2);
        if (bases_b[j + (upper << 2) + 3] < bases_a[i]) break;        
        int lower = (r >> 1);
        while (lower < upper) {
            int mid = (lower + upper) >> 1;
            if (bases_b[j + (mid << 2) + 3] >= bases_a[i]) upper = mid;
            else lower = mid + 1;
        }
        j += (lower << 2);

        __m128i bv_a = _mm_set_epi32(bases_a[i], bases_a[i], bases_a[i], bases_a[i]);
        __m128i bv_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i cmp_mask = _mm_cmpeq_epi32(bv_a, bv_b);
        int mask = _mm_movemask_ps((__m128)cmp_mask);
        if (mask != 0) {
            int p = __builtin_ctz(mask);
            states_c[size_c] = states_a[i] & states_b[j + p];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
        }
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int* prepare_byte_check_mask_dict2()
{
    int * mask = new int[65536];

    auto trans_c_s = [](const int c) -> int {
        switch (c) {
            case 0: return -1; // no match
            case 1: return 0;
            case 2: return 1;
            case 4: return 2;
            case 8: return 3;
            default: return 4; // multiple matches.
        }
    };

    for (int x = 0; x < 65536; ++x) {        
        int c0 = (x & 0xf), c1 = ((x >> 4) & 0xf);
        int c2 = ((x >> 8) & 0xf), c3 = ((x >> 12) & 0xf);
        int s0 = trans_c_s(c0), s1= trans_c_s(c1);
        int s2 = trans_c_s(c2), s3 = trans_c_s(c3);
        
        bool is_multiple_match = (s0 == 4) || (s1 == 4) ||
                (s2 == 4) || (s3 == 4);
        if (is_multiple_match) {
            mask[x] = -1;
            continue;
        }
        bool is_no_match = (s0 == -1) && (s1 == -1) &&
                (s2 == -1) && (s3 == -1);
        if (is_no_match) {
            mask[x] = -2;
            continue;
        }
        if (s0 == -1) s0 = 0; if (s1 == -1) s1 = 1;
        if (s2 == -1) s2 = 2; if (s3 == -1) s3 = 3;
        mask[x] = (s0) | (s1 << 2) | (s2 << 4) | (s3 << 6);        
    }

    return mask;
}
static const int *byte_check_mask_dict = prepare_byte_check_mask_dict2();

uint8_t * prepare_match_shuffle_dict2()
{
    uint8_t * dict = new uint8_t[4096];

    for (int x = 0; x < 256; ++x) {
        for (int i = 0; i < 4; ++i) {
            uint8_t c = (x >> (i << 1)) & 3; // c = 0, 1, 2, 3
            int pos = x * 16 + i * 4;
            for (uint8_t j = 0; j < 4; ++j)
                dict[pos + j] = c * 4 + j;
        }
    }

    return dict;
}
static const __m128i *match_shuffle_dict = (__m128i*)prepare_match_shuffle_dict2();

static const uint8_t byte_check_group_a_pi8[64] = {
    0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8, 12, 12, 12, 12,
    1, 1, 1, 1, 5, 5, 5, 5, 9, 9, 9, 9, 13, 13, 13, 13,
    2, 2, 2, 2, 6, 6, 6, 6, 10, 10, 10, 10, 14, 14, 14, 14,
    3, 3, 3, 3, 7, 7, 7, 7, 11, 11, 11, 11, 15, 15, 15, 15,
};
static const uint8_t byte_check_group_b_pi8[64] = {
    0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12,
    1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13,
    2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14,
    3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15,
};
static const __m128i *byte_check_group_a_order = (__m128i*)(byte_check_group_a_pi8);
static const __m128i *byte_check_group_b_order = (__m128i*)(byte_check_group_b_pi8);

int intersect_qfilter_uint_b4(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    // for (int x = 0; x < 65536; x += 16)
    //     _mm_prefetch((char*) (byte_check_mask_dict + x), _MM_HINT_T2);
    // for (int x = 0; x < 256; x += 4)
    //     _mm_prefetch((char*) (match_shuffle_dict + x), _MM_HINT_T0);
    // _mm_prefetch((char*) (byte_check_group_a_order), _MM_HINT_T0);
    // _mm_prefetch((char*) (byte_check_group_b_order), _MM_HINT_T0);
    
    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));
        
        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        }
      
        __m128i byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[0]);
        __m128i byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[0]);
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];
        if (__builtin_expect(ms_order == -1, 0)) {
            byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[1]);
            byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[1]);
            byte_check_mask = _mm_and_si128(byte_check_mask,
                    _mm_cmpeq_epi8(byte_group_a, byte_group_b));
            bc_mask = _mm_movemask_epi8(byte_check_mask);
            ms_order = byte_check_mask_dict[bc_mask];
            
            if (__builtin_expect(ms_order == -1, 0)) {
                byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[2]);
                byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[2]);
                byte_check_mask = _mm_and_si128(byte_check_mask,
                        _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                bc_mask = _mm_movemask_epi8(byte_check_mask);
                ms_order = byte_check_mask_dict[bc_mask];
                
                if (__builtin_expect(ms_order == -1, 0)) {
                    byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[3]);
                    byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[3]);
                    byte_check_mask = _mm_and_si128(byte_check_mask,
                            _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                    bc_mask = _mm_movemask_epi8(byte_check_mask);
                    ms_order = byte_check_mask_dict[bc_mask];                    
                }
            }
        }
        if (ms_order == -2) continue; // "no match" in this two block.

        __m128i sf_v_b = _mm_shuffle_epi8(v_b, match_shuffle_dict[ms_order]);
        __m128i cmp_mask = _mm_cmpeq_epi32(v_a, sf_v_b);

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), p);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_qfilter_uint_b4_v2(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);
 
    __m128i v_a = MM_LOAD_SI_128((__m128i*)set_a);
    __m128i v_b = MM_LOAD_SI_128((__m128i*)set_b);
    __m128i byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[0]);
    __m128i byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[0]);
    __m128i cmp_mask;
    while (i < qs_a && j < qs_b) {
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];

        if (__builtin_expect(ms_order != -2, 0)) {
            if (ms_order > 0) {
                 __m128i sf_v_b = _mm_shuffle_epi8(v_b, match_shuffle_dict[ms_order]);
                cmp_mask = _mm_cmpeq_epi32(v_a, sf_v_b);       
            } else {
                __m128i cmp_mask0 = _mm_cmpeq_epi32(v_a, v_b); // pairwise comparison
                __m128i rot1 = _mm_shuffle_epi32(v_b, cyclic_shift1);   // shuffling
                __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, rot1);
                __m128i rot2 = _mm_shuffle_epi32(v_b, cyclic_shift2);
                __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, rot2);
                __m128i rot3 = _mm_shuffle_epi32(v_b, cyclic_shift3);
                __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, rot3);
                cmp_mask = _mm_or_si128(
                        _mm_or_si128(cmp_mask0, cmp_mask1),
                        _mm_or_si128(cmp_mask2, cmp_mask3));
            } 

            int mask = _mm_movemask_ps((__m128)cmp_mask);
            __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
            _mm_storeu_si128((__m128i*)(set_c + size_c), p);

            size_c += _mm_popcnt_u32(mask);           
        }
    
        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        if (a_max <= b_max) {
            i += 4;
            v_a = MM_LOAD_SI_128((__m128i*)(set_a + i));
            byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[0]);
            _mm_prefetch((char*) (set_a + i + 16), _MM_HINT_T0);
        }
        if (a_max >= b_max) {
            j += 4;
            v_b = MM_LOAD_SI_128((__m128i*)(set_b + j));
            byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[0]);
            _mm_prefetch((char*) (set_b + j + 16), _MM_HINT_T0);
        }
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_qfilter_bsr_b4(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);
    
    while (i < qs_a && j < qs_b) {
        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i base_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));

        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        }       

        __m128i byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
        __m128i byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];
        if (__builtin_expect(ms_order == -1, 0)) {
            byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[1]);
            byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[1]);
            byte_check_mask = _mm_and_si128(byte_check_mask,
                    _mm_cmpeq_epi8(byte_group_a, byte_group_b));
            bc_mask = _mm_movemask_epi8(byte_check_mask);
            ms_order = byte_check_mask_dict[bc_mask];
            if (__builtin_expect(ms_order == -1, 0)) {
                byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[2]);
                byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[2]);
                byte_check_mask = _mm_and_si128(byte_check_mask,
                        _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                bc_mask = _mm_movemask_epi8(byte_check_mask);
                ms_order = byte_check_mask_dict[bc_mask];
                if (__builtin_expect(ms_order == -1, 0)) {
                    byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[3]);
                    byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[3]);
                    byte_check_mask = _mm_and_si128(byte_check_mask,
                            _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                    bc_mask = _mm_movemask_epi8(byte_check_mask);
                    ms_order = byte_check_mask_dict[bc_mask];
                }
            }
        }
        
        if (ms_order == -2) continue; // "no match" in this two block.

        __m128i sf_base_b = _mm_shuffle_epi8(base_b, match_shuffle_dict[ms_order]);
        __m128i sf_state_b = _mm_shuffle_epi8(state_b, match_shuffle_dict[ms_order]);
        __m128i cmp_mask = _mm_cmpeq_epi32(base_a, sf_base_b);        
        __m128i and_state = _mm_and_si128(state_a, sf_state_b);
        __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
        cmp_mask = _mm_andnot_si128(state_mask, cmp_mask);

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
        __m128i res_s = _mm_shuffle_epi8(and_state, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
        _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;      
}

int intersect_qfilter_bsr_b4_v2(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    __m128i base_a = _mm_lddqu_si128((__m128i*)bases_a + i);
    __m128i base_b = _mm_lddqu_si128((__m128i*)bases_b + j);
    __m128i byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
    __m128i byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
    __m128i cmp_mask, and_state;

    while (i < qs_a && j < qs_b) {
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];

        if (__builtin_expect(ms_order != -2, 0)) {
            __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
            __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));    
            if (ms_order > 0) {                
                __m128i sf_base_b = _mm_shuffle_epi8(base_b, match_shuffle_dict[ms_order]);
                __m128i sf_state_b = _mm_shuffle_epi8(state_b, match_shuffle_dict[ms_order]);
                cmp_mask = _mm_cmpeq_epi32(base_a, sf_base_b);                
                and_state = _mm_and_si128(state_a, sf_state_b);
                __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
                cmp_mask = _mm_andnot_si128(state_mask, cmp_mask);
            } else {
                __m128i cmp_mask0 = _mm_cmpeq_epi32(base_a, base_b);
                __m128i state_c0 = _mm_and_si128(
                        _mm_and_si128(state_a, state_b), cmp_mask0);
                __m128i base_sf1 = _mm_shuffle_epi32(base_b, cyclic_shift1);
                __m128i state_sf1 = _mm_shuffle_epi32(state_b, cyclic_shift1);
                __m128i cmp_mask1 = _mm_cmpeq_epi32(base_a, base_sf1);
                __m128i state_c1 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf1), cmp_mask1);
                __m128i base_sf2 = _mm_shuffle_epi32(base_b, cyclic_shift2);
                __m128i state_sf2 = _mm_shuffle_epi32(state_b, cyclic_shift2);
                __m128i cmp_mask2 = _mm_cmpeq_epi32(base_a, base_sf2);
                __m128i state_c2 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf2), cmp_mask2);
                __m128i base_sf3 = _mm_shuffle_epi32(base_b, cyclic_shift3);
                __m128i state_sf3 = _mm_shuffle_epi32(state_b, cyclic_shift3);
                __m128i cmp_mask3 = _mm_cmpeq_epi32(base_a, base_sf3);
                __m128i state_c3 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf3), cmp_mask3);
                and_state = _mm_or_si128(
                        _mm_or_si128(state_c0, state_c1),
                        _mm_or_si128(state_c2, state_c3)
                        );                
                __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
                cmp_mask = _mm_andnot_si128(state_mask, all_one_si128);
            } 

            int mask = _mm_movemask_ps((__m128)cmp_mask);
            __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
            __m128i res_s = _mm_shuffle_epi8(and_state, shuffle_mask[mask]);
            _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
            _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);

            size_c += _mm_popcnt_u32(mask);        
        }

        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        if (a_max <= b_max) {
            i += 4;
            base_a = MM_LOAD_SI_128((__m128i*)(bases_a + i));
            byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
            _mm_prefetch((char*) (bases_a + i + 16), _MM_HINT_T0);
            _mm_prefetch((char*) (states_a + i + 16), _MM_HINT_T0);
        }
        if (a_max >= b_max) {
            j += 4;
            base_b = MM_LOAD_SI_128((__m128i*)(bases_b + j));
            byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
            _mm_prefetch((char*) (bases_b + j + 16), _MM_HINT_T0);
            _mm_prefetch((char*) (states_b + j + 16), _MM_HINT_T0);
        }
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_shuffle_uint_b4(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));

        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        }

        __m128i cmp_mask0 = _mm_cmpeq_epi32(v_a, v_b); // pairwise comparison
        __m128i rot1 = _mm_shuffle_epi32(v_b, cyclic_shift1);   // shuffling
        __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, rot1);
        __m128i rot2 = _mm_shuffle_epi32(v_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, rot2);
        __m128i rot3 = _mm_shuffle_epi32(v_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, rot3);
        __m128i cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask0, cmp_mask1),
                _mm_or_si128(cmp_mask2, cmp_mask3));

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), p);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_shuffle_uint_b8(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 7);
    int qs_b = size_b - (size_b & 7);

    while (i < qs_a && j < qs_b) {
        __m128i v_a0 = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_a1 = _mm_lddqu_si128((__m128i*)(set_a + i + 4));
        __m128i v_b0 = _mm_lddqu_si128((__m128i*)(set_b + j));
        __m128i v_b1 = _mm_lddqu_si128((__m128i*)(set_b + j + 4));

        int a_max = set_a[i + 7];
        int b_max = set_b[j + 7];
        if (a_max == b_max) {
            i += 8;
            j += 8;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 8;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 8;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        }

        // a0 -- b0:
        __m128i cmp_mask0 = _mm_cmpeq_epi32(v_a0, v_b0); 
        __m128i rot1 = _mm_shuffle_epi32(v_b0, cyclic_shift1);   
        __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a0, rot1);
        __m128i rot2 = _mm_shuffle_epi32(v_b0, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a0, rot2);
        __m128i rot3 = _mm_shuffle_epi32(v_b0, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a0, rot3);

        // a0 -- b1:
        __m128i cmp_mask4 = _mm_cmpeq_epi32(v_a0, v_b1); 
        __m128i rot5 = _mm_shuffle_epi32(v_b1, cyclic_shift1);   
        __m128i cmp_mask5 = _mm_cmpeq_epi32(v_a0, rot5);
        __m128i rot6 = _mm_shuffle_epi32(v_b1, cyclic_shift2);
        __m128i cmp_mask6 = _mm_cmpeq_epi32(v_a0, rot6);
        __m128i rot7 = _mm_shuffle_epi32(v_b1, cyclic_shift3);
        __m128i cmp_mask7 = _mm_cmpeq_epi32(v_a0, rot7);

        __m128i cmp_maskx = _mm_or_si128(
                _mm_or_si128(
                    _mm_or_si128(cmp_mask0, cmp_mask1),
                    _mm_or_si128(cmp_mask2, cmp_mask3)
                ),
                _mm_or_si128(
                    _mm_or_si128(cmp_mask4, cmp_mask5),
                    _mm_or_si128(cmp_mask6, cmp_mask7)
                )
            );


        int maskx = _mm_movemask_ps((__m128)cmp_maskx);
        __m128i px = _mm_shuffle_epi8(v_a0, shuffle_mask[maskx]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), px);
        size_c += _mm_popcnt_u32(maskx);

        // a1 -- b0:
        cmp_mask0 = _mm_cmpeq_epi32(v_a1, v_b0); 
        rot1 = _mm_shuffle_epi32(v_b0, cyclic_shift1);   
        cmp_mask1 = _mm_cmpeq_epi32(v_a1, rot1);
        rot2 = _mm_shuffle_epi32(v_b0, cyclic_shift2);
        cmp_mask2 = _mm_cmpeq_epi32(v_a1, rot2);
        rot3 = _mm_shuffle_epi32(v_b0, cyclic_shift3);
        cmp_mask3 = _mm_cmpeq_epi32(v_a1, rot3);

        // a1 -- b1:
        cmp_mask4 = _mm_cmpeq_epi32(v_a1, v_b1); 
        rot5 = _mm_shuffle_epi32(v_b1, cyclic_shift1);   
        cmp_mask5 = _mm_cmpeq_epi32(v_a1, rot5);
        rot6 = _mm_shuffle_epi32(v_b1, cyclic_shift2);
        cmp_mask6 = _mm_cmpeq_epi32(v_a1, rot6);
        rot7 = _mm_shuffle_epi32(v_b1, cyclic_shift3);
        cmp_mask7 = _mm_cmpeq_epi32(v_a1, rot7);
        
        __m128i cmp_masky = _mm_or_si128(
                _mm_or_si128(
                    _mm_or_si128(cmp_mask0, cmp_mask1),
                    _mm_or_si128(cmp_mask2, cmp_mask3)
                ),
                _mm_or_si128(
                    _mm_or_si128(cmp_mask4, cmp_mask5),
                    _mm_or_si128(cmp_mask6, cmp_mask7)
                )
            );


        int masky = _mm_movemask_ps((__m128)cmp_masky);
        __m128i py = _mm_shuffle_epi8(v_a1, shuffle_mask[masky]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), py);
        size_c += _mm_popcnt_u32(masky);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c; 
}

int intersect_shuffle_bsr_b4(int* bases_a, int* states_a, int size_a,
            int* bases_b, int* states_b, int size_b,
            int* bases_c, int* states_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);
    
    while (i < qs_a && j < qs_b) {
        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i base_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));
        
        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        }        

        // shift0:
        __m128i cmp_mask0 = _mm_cmpeq_epi32(base_a, base_b);
        __m128i state_c0 = _mm_and_si128(
                _mm_and_si128(state_a, state_b), cmp_mask0);

        // shift1:
        __m128i base_sf1 = _mm_shuffle_epi32(base_b, cyclic_shift1);
        __m128i state_sf1 = _mm_shuffle_epi32(state_b, cyclic_shift1);
        __m128i cmp_mask1 = _mm_cmpeq_epi32(base_a, base_sf1);
        __m128i state_c1 = _mm_and_si128(
                _mm_and_si128(state_a, state_sf1), cmp_mask1);

        // shift2:
        __m128i base_sf2 = _mm_shuffle_epi32(base_b, cyclic_shift2);
        __m128i state_sf2 = _mm_shuffle_epi32(state_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(base_a, base_sf2);
        __m128i state_c2 = _mm_and_si128(
                _mm_and_si128(state_a, state_sf2), cmp_mask2);

        // shift3:
        __m128i base_sf3 = _mm_shuffle_epi32(base_b, cyclic_shift3);
        __m128i state_sf3 = _mm_shuffle_epi32(state_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(base_a, base_sf3);
        __m128i state_c3 = _mm_and_si128(
                _mm_and_si128(state_a, state_sf3), cmp_mask3);

        __m128i state_all = _mm_or_si128(
                _mm_or_si128(state_c0, state_c1),
                _mm_or_si128(state_c2, state_c3)
                );
        __m128i cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask0, cmp_mask1),
                _mm_or_si128(cmp_mask2, cmp_mask3)
                );
        __m128i state_mask = _mm_cmpeq_epi32(state_all, all_zero_si128);
        int mask = (_mm_movemask_ps((__m128)cmp_mask) &
                ~(_mm_movemask_ps((__m128)state_mask)));

        __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
        __m128i res_s = _mm_shuffle_epi8(state_all, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
        _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            bases_c[size_c] = bases_a[i];
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;    
}

constexpr int word_check_shuffle_a01 = _MM_SHUFFLE(1,1,0,0); 
constexpr int word_check_shuffle_a23 = _MM_SHUFFLE(3,3,2,2); 
constexpr int word_check_shuffle_b01 = _MM_SHUFFLE(1,0,1,0); 
constexpr int word_check_shuffle_b23 = _MM_SHUFFLE(3,2,3,2); 

int intersect_bmiss_uint_b4(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));

        // byte-wise check:
        __m128i byte_group_a0 = _mm_shuffle_epi8(v_a, byte_check_group_a_order[0]);
        __m128i byte_group_b0 = _mm_shuffle_epi8(v_b, byte_check_group_b_order[0]);
        __m128i byte_check_mask0 = _mm_cmpeq_epi8(byte_group_a0, byte_group_b0);
        int bc_mask0 = _mm_movemask_epi8(byte_check_mask0);
        __m128i byte_group_a1 = _mm_shuffle_epi8(v_a, byte_check_group_a_order[1]);
        __m128i byte_group_b1 = _mm_shuffle_epi8(v_b, byte_check_group_b_order[1]);
        __m128i byte_check_mask1 = _mm_cmpeq_epi8(byte_group_a1, byte_group_b1);
        int bc_mask1 = _mm_movemask_epi8(byte_check_mask1);
        if ((bc_mask0 & bc_mask1) == 0) {
            if (set_a[i + 3] == set_b[j + 3]) {
                i += 4;
                j += 4;
            } else if (set_a[i + 3] < set_b[j + 3]) {
                i += 4;
            } else {
                j += 4;
            }
            continue;
        }

        // word-wise check:
        __m128i word_group_a01 = _mm_shuffle_epi32(v_a, word_check_shuffle_a01);
        __m128i word_group_b01 = _mm_shuffle_epi32(v_b, word_check_shuffle_b01);
        __m128i word_group_a23 = _mm_shuffle_epi32(v_a, word_check_shuffle_a23);
        __m128i word_group_b23 = _mm_shuffle_epi32(v_b, word_check_shuffle_b23);

        __m128i work_check_mask0101 = _mm_cmpeq_epi32(word_group_a01, word_group_b01);
        __m128i work_check_mask0123 = _mm_cmpeq_epi32(word_group_a01, word_group_b23);
        __m128i work_check_mask01 = _mm_or_si128(work_check_mask0101, work_check_mask0123);

        __m128i work_check_mask2301 = _mm_cmpeq_epi32(word_group_a23, word_group_b01);
        __m128i work_check_mask2323 = _mm_cmpeq_epi32(word_group_a23, word_group_b23);
        __m128i work_check_mask23 = _mm_or_si128(work_check_mask2301, work_check_mask2323);

        int wc_mask01 = _mm_movemask_ps((__m128)work_check_mask01);
        if ((wc_mask01 & 0x3) != 0) set_c[size_c++] = set_a[i];
        if ((wc_mask01 & 0xc) != 0) set_c[size_c++] = set_a[i + 1];
        int wc_mask23 = _mm_movemask_ps((__m128)work_check_mask23);
        if ((wc_mask23 & 0x3) != 0) set_c[size_c++] = set_a[i + 2];
        if ((wc_mask23 & 0xc) != 0) set_c[size_c++] = set_a[i + 3];

        if (set_a[i + 3] == set_b[j + 3]) {
            i += 4;
            j += 4;
        } else if (set_a[i + 3] < set_b[j + 3]) {
            i += 4;
        } else {
            j += 4;
        }
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

static const uint8_t bmiss_sttni_bc_array[32] = {
    0, 1, 4, 5, 8, 9, 12, 13, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 0, 1, 4, 5, 8, 9, 12, 13,
};
static const __m128i *BMISS_BC_ORD = (__m128i*)(bmiss_sttni_bc_array);

int intersect_bmiss_uint_sttni_b8(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 7);
    int qs_b = size_b - (size_b & 7); 

    while (i < qs_a && j < qs_b) {
        __m128i v_a0 = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_a1 = _mm_lddqu_si128((__m128i*)(set_a + i + 4));
        __m128i v_b0 = _mm_lddqu_si128((__m128i*)(set_b + j));
        __m128i v_b1 = _mm_lddqu_si128((__m128i*)(set_b + j + 4));

        // byte-wise check by STTNI:
        __m128i byte_group_a0 = _mm_shuffle_epi8(v_a0, BMISS_BC_ORD[0]);
        __m128i byte_group_a1 = _mm_shuffle_epi8(v_a1, BMISS_BC_ORD[1]);
        __m128i byte_group_a = _mm_or_si128(byte_group_a0, byte_group_a1);
        __m128i byte_group_b0 = _mm_shuffle_epi8(v_b0, BMISS_BC_ORD[0]);
        __m128i byte_group_b1 = _mm_shuffle_epi8(v_b1, BMISS_BC_ORD[1]);
        __m128i byte_group_b = _mm_or_si128(byte_group_b0, byte_group_b1);

        __m128i bc_mask = _mm_cmpestrm(byte_group_b, 8, byte_group_a, 8,
                _SIDD_UWORD_OPS|_SIDD_CMP_EQUAL_ANY|_SIDD_BIT_MASK);
        int r = _mm_extract_epi32(bc_mask, 0);

        // word-wise check:
        while (r) {
            int p = _mm_popcnt_u32((~r) & (r - 1));
            r &= (r - 1);
            __m128i wc_a = _mm_set_epi32(set_a[i + p], set_a[i + p],
                            set_a[i + p], set_a[i + p]);
            if (!_mm_test_all_zeros(_mm_cmpeq_epi32(wc_a, v_b0),
                                all_one_si128) ||
                !_mm_test_all_zeros(_mm_cmpeq_epi32(wc_a, v_b1),
                                all_one_si128))
                set_c[size_c++] = set_a[i + p];
        }

        if (set_a[i + 7] == set_b[j + 7]) {
            i += 8;
            j += 8;
        } else if (set_a[i + 7] < set_b[j + 7]) {
            i += 8;
        } else {
            j += 8;
        }
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c; 
}

uint8_t* prepare_shuffling_dict_u16()
{
    uint8_t* mask = new uint8_t[4096];
    memset(mask, 0xff, sizeof(uint8_t) * 4096);
    int size = 0;
    for (int i = 0; i < 256; ++i) {
        int counter = 0;
        for (int j = 0; j < 8; ++j) {
            if (i & (1 << j)) {
                mask[size + counter * 2    ] = 2 * j;
                mask[size + counter * 2 + 1] = 2 * j + 1;
                ++counter;             
             }              
        }
        size += 16;
    }
    return mask;
}
static const __m128i *shuffle_mask16 = (__m128i *)prepare_shuffling_dict_u16();

int intersect_hierainter_uint_sttni(int *set_a, int size_a,
        int *set_b, int size_b, int *set_c)
{
    uint16_t *hi_set_a = NULL, *hi_set_b = NULL, *hi_set_c = NULL;
    align_malloc((void**)&hi_set_a, 32, sizeof(uint16_t) * size_a * 3);
    align_malloc((void**)&hi_set_b, 32, sizeof(uint16_t) * size_b * 3);
    align_malloc((void**)&hi_set_c, 32, sizeof(uint16_t) * std::min(size_a, size_b) * 3);
    int hi_size_a = hierainter_offline_partition(set_a, size_a, hi_set_a);
    int hi_size_b = hierainter_offline_partition(set_b, size_b, hi_set_b);

//    struct timeval time_start;
//    struct timeval time_end;
//    gettimeofday(&time_start, NULL);

    int hi_size_c = hierainter_online_intersect_high16bit(hi_set_a, hi_size_a,
                hi_set_b, hi_size_b, hi_set_c);

//    gettimeofday(&time_end, NULL);
//    double run_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 +
//                (time_end.tv_usec - time_start.tv_usec) / 1000.0;
//    if (run_time > 10000.0) printf("hierainter_online_time=%.3fs\n", run_time / 1000.0);
//    else printf("hierainter_online_time=%.3fms\n", run_time);

    int size_c = hierainter_offline_combine(hi_set_c, hi_size_c, set_c);
    return size_c;
}

int hierainter_offline_partition(int *set_a, int size_a, uint16_t *hi_set)
{
    uint16_t high = 0;
    int partition_len = 0, partition_pos = 1;
    int hi_size = 0;

    for (int i = 0; i < size_a; ++i) {
        uint16_t chigh = (uint16_t)(set_a[i] >> 16); // upper 16-bit part
        uint16_t clow = (uint16_t)(set_a[i] & 0xffff); // lower 16-bit part
        if (chigh == high && i != 0) { // add element to the current partition
            hi_set[hi_size++] = clow;
            partition_len++;
        } else { // start new partition
            hi_set[hi_size++] = chigh; // partition prefix
            hi_set[hi_size++] = 0;     // reserve place for partition size
            hi_set[hi_size++] = clow;  // write the first element
            hi_set[partition_pos] = partition_len;
            partition_len = 1; // reset counters
            partition_pos = hi_size - 2;
            high = chigh;            
        }
    }
    hi_set[partition_pos] = partition_len;

    return hi_size;
}

int hierainter_offline_combine(uint16_t *hi_set, int hi_size, int *set_a)
{
    int i = 0, size_a = 0;
    while (i < hi_size) {
        int chigh = (int)hi_set[i++]; chigh <<= 16;
        int par_len = (int)hi_set[i++];
        for (int j = 0; j < par_len; ++j) {
            int clow = (int)hi_set[i + j];
            set_a[size_a++] = (chigh | clow);
        }
        i += par_len;
    }
    return size_a;
}

int hierainter_online_intersect_low16bit(uint16_t *set_a, int size_a,
        uint16_t *set_b, int size_b, uint16_t *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 7);
    int qs_b = size_b - (size_b & 7);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));
 
        __m128i res_v = _mm_cmpestrm(v_b, 8, v_a, 8,
                _SIDD_UWORD_OPS|_SIDD_CMP_EQUAL_ANY|_SIDD_BIT_MASK);
        int r = _mm_extract_epi32(res_v, 0);
        __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask16[r]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), p);
        size_c += _mm_popcnt_u32(r);
 
        uint16_t a_max = _mm_extract_epi16(v_a, 7);
        uint16_t b_max = _mm_extract_epi16(v_b, 7);
        i += (a_max <= b_max) * 8;
        j += (a_max >= b_max) * 8;
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int hierainter_online_intersect_high16bit(uint16_t *par_a, int size_a,
        uint16_t *par_b, int size_b, uint16_t *par_c)
{
    int i = 0, j = 0, size_c = 0;

    while (i < size_a && j < size_b) {
        if (par_a[i] < par_b[j]) {
            i += par_a[i + 1] + 2;
        } else if (par_b[j] < par_a[i]) {
            j += par_b[j + 1] + 2;
        } else {
            par_c[size_c++] = par_a[i];  // write partition prefix            
            uint16_t partition_len = hierainter_online_intersect_low16bit(
                    par_a + i + 2, par_a[i + 1], par_b + j + 2, par_b[j + 1],
                    par_c + size_c + 1);
            par_c[size_c++] = partition_len;
            size_c += partition_len;
            i += par_a[i + 1] + 2;
            j += par_b[j + 1] + 2;
        }
    }
    return size_c;
}


const int BSR_WIDTH = sizeof(int) * 8;
const int BSR_SHIFT = __builtin_ctzll(BSR_WIDTH);
const int BSR_MASK = BSR_WIDTH - 1;
int offline_uint_trans_bsr(int *set_a, int size_a, int *bases_a, int *states_a)
{
    int cnt = -1;
    for (int i = 0; i < size_a; ++i) {
        int u = set_a[i];
        int u_base = (u >> BSR_SHIFT);
        int u_bit = (1 << (u & BSR_MASK));
        if (cnt == -1 || bases_a[cnt] != u_base) {
            bases_a[++cnt] = u_base;
            states_a[cnt] = u_bit;
        } else {
            states_a[cnt] |= u_bit;
        }
    }
    return ++cnt;    
}

int offline_bsr_trans_uint(int *bases_a, int *states_a, int size_a, int *set_a)
{
    int cnt = 0;
    for (int i = 0; i < size_a; ++i) {
        int u_high = (bases_a[i] << BSR_SHIFT);
        int state = states_a[i];
        while (state) {
            int u = (u_high | __builtin_ctz(state));
            set_a[cnt++] = u;
            state &= (state - 1);
        }
    }
    return cnt;
}