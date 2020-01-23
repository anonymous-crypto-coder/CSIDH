#ifndef FP_H
#define FP_H

#include "u512.h"

/* fp is in the Montgomery domain, so interpreting that
   as an integer should never make sense.
   enable compiler warnings when mixing up u512 and fp. */
typedef struct fp {
    u512 x;
} fp;

extern const fp fp_0;
extern const fp fp_1;

void fp_set(fp *x, uint64_t y);         /* set value of x to y */
void fp_cswap(fp *x, fp *y, bool c);    /* conditional swap */

void fp_enc(fp *x, u512 const *y); /* encode to Montgomery representation */
void fp_dec(u512 *x, fp const *y); /* decode from Montgomery representation */

/* a '2' function performs the indicated operation on the 2 inputs and stores the result in the first */
void fp_add2(fp *x, fp const *y);
void fp_sub2(fp *x, fp const *y);
void fp_mul2(fp *x, fp const *y);

/* a '3' function performs the indicated  operation on the last 2 inputs and stores the result in the first*/
void fp_add3(fp *x, fp const *y, fp const *z);
void fp_sub3(fp *x, fp const *y, fp const *z);
void fp_mul3(fp *x, fp const *y, fp const *z);

void fp_sq1(fp *x);                 /* squares x, updates x to result */
void fp_sq2(fp *x, fp const *y);    /* squares y, updates x to result */
void fp_inv(fp *x);                 /* inverts x, updates x to result */
bool fp_issquare(fp const *x);      /* checks if x is a square */

void fp_random(fp *x);



#endif
