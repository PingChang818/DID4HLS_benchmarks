/*
Implementation based on algorithm described in:
A. Danalis, G. Marin, C. McCurdy, J. S. Meredith, P. C. Roth, K. Spafford, V. Tipparaju, and J. S. Vetter.
The scalable heterogeneous computing (shoc) benchmark suite.
In Proceedings of the 3rd Workshop on General-Purpose Computation on Graphics Processing Units, 2010
*/
// arr_nest = [-1, -1, 1, -1, 3, -1, -1, 6, -1, 8]
// arr_unroll = [32, 16, 3, 8, 14, 32, 8, 15, 8, 3]

#include "sort.h"

void local_scan (int bucket[BUCKETSIZE]) {
    int radixID, i, bucket_indx;
    loop_3: for (radixID = 0; radixID < SCAN_RADIX; radixID++) {
        loop_4: for (i = 1; i < SCAN_BLOCK; i++){
            bucket_indx = radixID * SCAN_BLOCK + i;
            bucket[bucket_indx] += bucket[bucket_indx - 1];
        }
    }
}

void sum_scan (int sum[SCAN_RADIX], int bucket[BUCKETSIZE]) {
    int radixID, bucket_indx;
    sum[0] = 0;
    loop_5: for (radixID = 1; radixID < SCAN_RADIX; radixID++) {
        bucket_indx = radixID * SCAN_BLOCK - 1;
        sum[radixID] = sum[radixID - 1] + bucket[bucket_indx];
    }
}

void last_step_scan (int bucket[BUCKETSIZE], int sum[SCAN_RADIX]) {
    int radixID, i, bucket_indx;
    loop_6: for (radixID = 0; radixID < SCAN_RADIX; radixID++) {
        loop_7: for (i = 0; i < SCAN_BLOCK; i++) {
            bucket_indx = radixID * SCAN_BLOCK + i;
            bucket[bucket_indx] = bucket[bucket_indx] + sum[radixID];
         }
    }
}

void init (int bucket[BUCKETSIZE]) {
    int i;
    loop_0: for (i = 0; i < BUCKETSIZE; i++) {
        bucket[i] = 0;
    }
}

void hist (int bucket[BUCKETSIZE], int a[SIZE], int exp) {
    int blockID, i, bucket_indx, a_indx;
    blockID = 0;
    loop_1: for (blockID = 0; blockID < NUMOFBLOCKS; blockID++) {
        loop_2: for (i = 0; i < 4; i++) {
            a_indx = blockID * ELEMENTSPERBLOCK + i;
            bucket_indx = ((a[a_indx] >> exp) & 0x3) * NUMOFBLOCKS + blockID + 1;
            bucket[bucket_indx]++;
        }
    }
}

void update (int b[SIZE], int bucket[BUCKETSIZE], int a[SIZE], int exp) {
    int i, blockID, bucket_indx, a_indx;
    blockID = 0;

    loop_8: for (blockID = 0; blockID < NUMOFBLOCKS; blockID++) {
        loop_9: for (i = 0; i < 4; i++) {
            bucket_indx = ((a[blockID * ELEMENTSPERBLOCK + i] >> exp) & 0x3) * NUMOFBLOCKS + blockID;
            a_indx = blockID * ELEMENTSPERBLOCK + i;
            b[bucket[bucket_indx]] = a[a_indx];
            bucket[bucket_indx]++;
        }
    }
}

void sort (int a[SIZE], int b[SIZE], int bucket[BUCKETSIZE], int sum[SCAN_RADIX]) {
    int exp=0;
    int valid_buffer=0;
    #define BUFFER_A 0
    #define BUFFER_B 1

    for (exp = 0; exp < 32; exp += 2) {
        init(bucket);
        if(valid_buffer == BUFFER_A){
            hist(bucket, a, exp);
        }else{
            hist(bucket, b, exp);
        }

        local_scan(bucket);
        sum_scan(sum, bucket);
        last_step_scan(bucket, sum);

        if(valid_buffer == BUFFER_A){
            update(b, bucket, a, exp);
            valid_buffer = BUFFER_B;
        }else{
            update(a, bucket, b, exp);
            valid_buffer = BUFFER_A;
        }
    }
    // If trip count is even, buffer A will be valid at the end.
}
