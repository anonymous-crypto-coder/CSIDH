#ifndef CSIDH_H
#define CSIDH_H

#include "u512.h"
#include "fp.h"
#include "mont.h"
#include "helper.h"

uint8_t elligator_index;
void assign_bounds(int8_t *bound);

void action(public_key *out, public_key const *in, private_key const *priv,
		int8_t const *max_exponent);
bool csidh(public_key *out, public_key const *in, private_key const *priv,
		int8_t const *max_exponent);

extern void fp_print(fp const *x);

//int e_temp[num_primes];
//int f_temp[num_primes];

#endif
