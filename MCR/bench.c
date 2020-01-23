#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "u512.h"
#include "fp.h"
#include "mont.h"
#include "csidh.h"
#include "cycle.h"

#include <inttypes.h>

static __inline__ uint64_t rdtsc(void) {
	uint32_t hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return lo | (uint64_t) hi << 32;
}

unsigned long its = 1024;



int main() {
	clock_t t0, t1, time = 0;
	uint64_t c0, c1, cycles = 0;
	uint64_t allticks = 0;
	ticks ticks1, ticks2;
	//int e_counter[its][num_primes];
	//int f_counter[its][num_primes];

	//uint8_t num_batches = 5;
	//uint8_t my = 11;

	int8_t exponent_bounds[num_primes] = { 8, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17, 16, 15, 13, 13, 13, 13, 13, 12, 12, 11, 11, 11, 10, 11, 10, 10, 10, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 6, 7, 7, 7, 7, 6, 7, 7, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 6, 5, 5, 5, 6, 4};

	//unsigned int num_isogenies = 763;

	// calculate inverses for "elligatoring"
	// create inverse of u^2 - 1 : from 2 - 11
	for (int i = 2; i <= 10; i++) {
		fp_set(&invs_[i - 2], i);
		fp_sq1(&invs_[i - 2]);
		fp_sub2(&invs_[i - 2], &fp_1);
		fp_inv(&invs_[i - 2]);
	}

	private_key priv;
	public_key pub = base;
	csidh_private(&priv, exponent_bounds);
	action(&pub, &pub, &priv, exponent_bounds);

	for (unsigned long i = 0; i < its; ++i) {

		csidh_private(&priv, exponent_bounds);

		t0 = clock();
		c0 = rdtsc();
		ticks1 = getticks();
		/**************************************/
		assert(validate(&pub));
		action(&pub, &pub, &priv, exponent_bounds);
		/**************************************/
		ticks2 = getticks();
		allticks = allticks + elapsed(ticks2, ticks1);
		c1 = rdtsc();
		t1 = clock();
		cycles += c1 - c0;
		time += t1 - t0;

		/*for(int j = 0; j < num_primes; j++){
			e_counter[i][j] = e_temp[j];
			f_counter[i][j] = f_temp[j];
		}*/
	}

	printf("\niterations: %lu\n", its);
	printf("clock cycles: %" PRIu64 " (rdtsc)\n", (uint64_t) cycles / its);
	printf("clock cycles: %" PRIu64 " (getticks)\n", (uint64_t) allticks / its);

	printf("wall-clock time: %.3lf ms\n", 1000. * time / CLOCKS_PER_SEC / its);

	/*
	float e_vec[num_primes], f_vec[num_primes];
	for (int i = 0; i < num_primes; i++){
		e_vec[i] = 0;
		f_vec[i] = 0;
		for (unsigned long j = 0; j < its; j++){
			e_vec[i] += e_counter[j][i];
			f_vec[i] += f_counter[j][i];
		}
		e_vec[i] = e_vec[i] / its;
		f_vec[i] = f_vec[i] / its;
	}
	printf("\ne vs f: ");
	for (int i = 0; i < num_primes; i++)
		printf("\n  %f  %f", e_vec[i], f_vec[i]);
	//printf("\nf vec: ");
	//for (int i = 0; i < num_primes; i++)
	//	printf("\n  %f", f_vec[i]);
	printf("\n\n");*/
}

