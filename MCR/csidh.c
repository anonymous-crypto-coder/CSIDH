#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "csidh.h"
#include "rng.h"


const unsigned primes[num_primes] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587};


void assign_bounds(int8_t *bound){
	int8_t vals[74] = {8, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17, 16, 15, 13, 13, 13, 13, 13, 12, 12, 11, 11, 11, 10, 11, 10, 10, 10, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 6, 7, 7, 7, 7, 6, 7, 7, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 6, 5, 5, 5, 6, 4};
	for (int i = 0; i < 74; i++)
		bound[i] = vals[i];
}


void strat11(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[30] = {72, 46, 45, 52, 49, 61, 62, 63, 65, 22, 67, 16, 20, 40, 64, 21, 23, 8, 19, 12, 11, 6, 58, 24, 25, 38, 42, 13, 39, 18};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0xf85189b5512c8e51, 0xeacc7a46830953c7, 0x6eba44e912a13564, 0x8868a6d9283903f7, 0x39fe450238b, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 30; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[4];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x2f49e9b189f57a43, 0x3415ef461a26b20c, 0x8ab9147, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x42c1614f6c1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xb3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1e24b0cf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xc1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2c9aeb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x6edd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x119, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(58, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 281) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 281, x);
		e[58] = ec - (1 ^ x);
		counter[58] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x125d67ed68f01e9, 0xd7c030d88a6, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x46958ab, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xc43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x66c7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x13d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(64, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 317) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 317, x);
		e[64] = ec - (1 ^ x);
		counter[64] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x753d4eb7c42dcce5, 0xbf59, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x198401, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x52af, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x15b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(67, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 347) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 347, x);
		e[67] = ec - (1 ^ x);
		counter[67] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xa9af8eefbf464727, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x14b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(65, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 331) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 331, x);
		e[65] = ec - (1 ^ x);
		counter[65] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x15c2a2b4bdb6b5f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(63, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 313) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 313, x);
		e[63] = ec - (1 ^ x);
		counter[63] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11e9799599b19, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(62, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 311) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 311, x);
		e[62] = ec - (1 ^ x);
		counter[62] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xeefb7f0a83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(61, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 307) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 307, x);
		e[61] = ec - (1 ^ x);
		counter[61] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x10692ac8b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(49, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 233) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 233, x);
		e[49] = ec - (1 ^ x);
		counter[49] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x10bcdb1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(52, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 251) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 251, x);
		e[52] = ec - (1 ^ x);
		counter[52] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x144eb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(45, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 211) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 211, x);
		e[45] = ec - (1 ^ x);
		counter[45] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x175, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(46, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 223) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 223, x);
		e[46] = ec - (1 ^ x);
		counter[46] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(72, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 373) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 373, x);
		e[72] = ec - (1 ^ x);
		counter[72] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat12(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[18] = {71, 48, 51, 47, 15, 66, 69, 29, 26, 60, 41, 59, 27, 43, 35, 68, 33, 32};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0xdebcf84a4f39043f, 0x8829a6b5d780375, 0x51fcd9a5ebf56675, 0xa8c1f44f446d38dc, 0x127a9c7d9167ffbf, 0xd167d04b3dddf, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 18; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[2];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xd76ef235af75423, 0x482614b19426a, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xcb21, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x15d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(68, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 349) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 349, x);
		e[68] = ec - (1 ^ x);
		counter[68] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xfd611f0f82fe43cd, 0x144cb9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x5cb9bb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x787f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x11b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(59, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 283) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 283, x);
		e[59] = ec - (1 ^ x);
		counter[59] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc5acfae2e5ca6db7, 0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x125, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(60, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 293) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 293, x);
		e[60] = ec - (1 ^ x);
		counter[60] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x5324a634997b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xaf1f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x161, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(69, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 353) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 353, x);
		e[69] = ec - (1 ^ x);
		counter[69] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x3f28c11b0b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(66, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 337) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 337, x);
		e[66] = ec - (1 ^ x);
		counter[66] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1350e9b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xe3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(47, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 227) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 227, x);
		e[47] = ec - (1 ^ x);
		counter[47] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1484b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(51, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 241) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 241, x);
		e[51] = ec - (1 ^ x);
		counter[51] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x16f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(48, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 229) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 229, x);
		e[48] = ec - (1 ^ x);
		counter[48] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(71, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 367) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 367, x);
		e[71] = ec - (1 ^ x);
		counter[71] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat13(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[26] = {73, 50, 55, 57, 53, 54, 36, 56, 44, 34, 14, 70, 28, 31, 37, 17, 30, 9, 7, 10, 4, 5, 3, 2, 0, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x3979b8bcde0a14f, 0xf4fe27a521743cd7, 0x5fd9632741e4928a, 0x1ba4365a3550e79f, 0xddd3aeb36ca8aa28, 0xda0d831, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 26; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[8];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x43e184ffbb277e03, 0x291deb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xe7a81d1bb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x16cdb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1e1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[7] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 8, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[7];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x54cdaf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xa7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x9e77, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x167, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(70, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 359) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 359, x);
		e[70] = ec - (1 ^ x);
		counter[70] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1992bbfdefc21d9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xd2a9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(56, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 271) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 271, x);
		e[56] = ec - (1 ^ x);
		counter[56] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x27185296855, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x107, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(54, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 263) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 263, x);
		e[54] = ec - (1 ^ x);
		counter[54] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x26f161355, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(53, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 257) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 257, x);
		e[53] = ec - (1 ^ x);
		counter[53] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x23fd941, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(57, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 277) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 277, x);
		e[57] = ec - (1 ^ x);
		counter[57] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x22405, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(55, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 269) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 269, x);
		e[55] = ec - (1 ^ x);
		counter[55] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x24b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(50, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 239) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 239, x);
		e[50] = ec - (1 ^ x);
		counter[50] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(73, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 587) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 587, x);
		e[73] = ec - (1 ^ x);
		counter[73] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat21(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[35] = {71, 49, 50, 51, 52, 45, 58, 69, 67, 26, 54, 22, 20, 30, 55, 15, 16, 2, 41, 19, 9, 11, 6, 56, 24, 25, 0, 1, 28, 43, 14, 31, 10, 18, 7};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x8ba77d850c54eecd, 0xdbda506135335582, 0x66a53f8234d4566, 0xbd283c99e9ee5074, 0x3202842, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 35; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[4];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xa14cc00ef1bcadf5, 0xb9f69071fcb7c035, 0x5d0435a3, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x2d609f04546b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x13cd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x11cd32193, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xc5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2854443, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x6aeb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(56, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 271) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 271, x);
		e[56] = ec - (1 ^ x);
		counter[56] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x3ae021299cc1298b, 0xb7953189f92, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4d26d76b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x8d7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x676855, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3dff, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(55, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 269) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 269, x);
		e[55] = ec - (1 ^ x);
		counter[55] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x9a3ddda4294c2f19, 0xcb6f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1c3741, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x5b6f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x107, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(54, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 263) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 263, x);
		e[54] = ec - (1 ^ x);
		counter[54] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x671525f494610191, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x15b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(67, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 347) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 347, x);
		e[67] = ec - (1 ^ x);
		counter[67] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x104694369937e31, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(69, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 353) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 353, x);
		e[69] = ec - (1 ^ x);
		counter[69] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xed3e30a910d9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(58, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 281) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 281, x);
		e[58] = ec - (1 ^ x);
		counter[58] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11fd6fabc23, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(45, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 211) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 211, x);
		e[45] = ec - (1 ^ x);
		counter[45] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x12592d8f9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(52, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 251) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 251, x);
		e[52] = ec - (1 ^ x);
		counter[52] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x137d889, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(51, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 241) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 241, x);
		e[51] = ec - (1 ^ x);
		counter[51] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x14e07, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(50, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 239) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 239, x);
		e[50] = ec - (1 ^ x);
		counter[50] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x16f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(49, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 233) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 233, x);
		e[49] = ec - (1 ^ x);
		counter[49] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(71, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 367) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 367, x);
		e[71] = ec - (1 ^ x);
		counter[71] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat22(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[21] = {72, 46, 47, 65, 68, 66, 62, 29, 34, 53, 27, 8, 33, 13, 5, 70, 23, 42, 44, 32, 17};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x8a694eccb4795b55, 0xadcf932213979c78, 0x37ac7484b359396d, 0x1cf6d6fecb91b866, 0x2fda25959928c959, 0x1b7c9690b1, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 21; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[4];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x857a91d4d4cadef3, 0x14e77e5f909a, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4fb7d231, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x8b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x668d47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x8807, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x167, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(70, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 359) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 359, x);
		e[70] = ec - (1 ^ x);
		counter[70] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xa74b8d67e1446461, 0xed, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xc6559, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x101, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x6d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(53, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 257) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 257, x);
		e[53] = ec - (1 ^ x);
		counter[53] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x29c88869a638f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x9a49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x137, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(62, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 311) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 311, x);
		e[62] = ec - (1 ^ x);
		counter[62] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1fbd8e5dedf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(66, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 337) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 337, x);
		e[66] = ec - (1 ^ x);
		counter[66] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x17484b16b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(68, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 349) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 349, x);
		e[68] = ec - (1 ^ x);
		counter[68] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1201c61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(65, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 331) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 331, x);
		e[65] = ec - (1 ^ x);
		counter[65] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x144eb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(47, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 227) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 227, x);
		e[47] = ec - (1 ^ x);
		counter[47] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x175, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(46, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 223) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 223, x);
		e[46] = ec - (1 ^ x);
		counter[46] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(72, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 373) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 373, x);
		e[72] = ec - (1 ^ x);
		counter[72] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat23(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[18] = {73, 48, 59, 60, 61, 21, 63, 39, 64, 40, 37, 57, 38, 35, 36, 12, 4, 3};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0xb1d8bb2ed0c20bb1, 0xd5cbced17e832b27, 0xfca6637a45311e8e, 0x5369a0be53219688, 0x8a42ce2675da7946, 0x78669f43508d87c, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 18; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[5];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xeb3c789bb0e569a5, 0x7e1f9f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x72cd0d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xa3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xbb31, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x115, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(57, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 277) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 277, x);
		e[57] = ec - (1 ^ x);
		counter[57] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xdcd4f462c9e53273, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xe021, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x13d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(64, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 317) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 317, x);
		e[64] = ec - (1 ^ x);
		counter[64] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1024fd37d7e49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x139, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(63, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 313) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 313, x);
		e[63] = ec - (1 ^ x);
		counter[63] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x2985dd4c1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x133, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(61, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 307) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 307, x);
		e[61] = ec - (1 ^ x);
		counter[61] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x244786d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(60, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 293) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 293, x);
		e[60] = ec - (1 ^ x);
		counter[60] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x20d17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(59, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 283) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 283, x);
		e[59] = ec - (1 ^ x);
		counter[59] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x24b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(48, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 229) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 229, x);
		e[48] = ec - (1 ^ x);
		counter[48] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(73, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 587) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 587, x);
		e[73] = ec - (1 ^ x);
		counter[73] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat31(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[25] = {71, 46, 52, 51, 50, 48, 70, 59, 62, 28, 63, 15, 24, 43, 58, 18, 19, 9, 42, 56, 20, 21, 35, 38, 37};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x6d4522c8097e5d83, 0x9f5431875b1b1da0, 0x59686343e748db7b, 0x9870e9fba4fb89c3, 0xf8f13f1fb095a5fe, 0x18, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 25; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[3];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x253a6a82c2e146d9, 0x9b936e4b87c93eb6, 0x238e, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xb3cbcf7fb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10a0e847, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1b1d33, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x53a1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(56, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 271) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 271, x);
		e[56] = ec - (1 ^ x);
		counter[56] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x9e0ceca605e77b61, 0x11869e5f6984, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x2b0ebb9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x4def, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x119, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(58, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 281) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 281, x);
		e[58] = ec - (1 ^ x);
		counter[58] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xe9200fd01c1df783, 0xccdb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1c75cf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x4823, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x139, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(63, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 313) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 313, x);
		e[63] = ec - (1 ^ x);
		counter[63] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x7e07650b4848fde5, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x137, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(62, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 311) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 311, x);
		e[62] = ec - (1 ^ x);
		counter[62] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x15994b5dcf9ecff, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(59, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 283) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 283, x);
		e[59] = ec - (1 ^ x);
		counter[59] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xf66e529f00a9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(70, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 359) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 359, x);
		e[70] = ec - (1 ^ x);
		counter[70] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1137c72b875, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(48, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 229) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 229, x);
		e[48] = ec - (1 ^ x);
		counter[48] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x12714d4db, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(50, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 239) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 239, x);
		e[50] = ec - (1 ^ x);
		counter[50] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x139728b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(51, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 241) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 241, x);
		e[51] = ec - (1 ^ x);
		counter[51] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x13fb1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(52, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 251) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 251, x);
		e[52] = ec - (1 ^ x);
		counter[52] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x16f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(46, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 223) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 223, x);
		e[46] = ec - (1 ^ x);
		counter[46] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(71, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 367) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 367, x);
		e[71] = ec - (1 ^ x);
		counter[71] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat32(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[25] = {72, 49, 45, 47, 67, 61, 60, 26, 57, 22, 36, 66, 14, 23, 39, 40, 65, 16, 25, 44, 41, 10, 30, 13, 17};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x33a3d24b7bb859b3, 0x1bd0027bc323faed, 0x3ce8f814ffd22f41, 0xf789e2d7d6330665, 0xc8e05a6c88ac809e, 0xa4, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 25; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[3];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x3c59406d6fa29887, 0x25bf298aa1e3ca0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x2a8f6773655, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x180d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x18aaeccf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xbf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1fbbb9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x4edf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x14b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(65, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 331) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 331, x);
		e[65] = ec - (1 ^ x);
		counter[65] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xb068d8c16e6f9215, 0x2e3624, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x127c105f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1a6fa5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x45c5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x151, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(66, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 337) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 337, x);
		e[66] = ec - (1 ^ x);
		counter[66] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc0ef88e91676a903, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x604d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x115, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(57, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 277) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 277, x);
		e[57] = ec - (1 ^ x);
		counter[57] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1934fdbe25cd5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x125, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(60, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 293) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 293, x);
		e[60] = ec - (1 ^ x);
		counter[60] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1504fef39d7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(61, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 307) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 307, x);
		e[61] = ec - (1 ^ x);
		counter[61] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xf81d7635, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(67, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 347) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 347, x);
		e[67] = ec - (1 ^ x);
		counter[67] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x117d007, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(47, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 227) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 227, x);
		e[47] = ec - (1 ^ x);
		counter[47] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1537d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(45, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 211) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 211, x);
		e[45] = ec - (1 ^ x);
		counter[45] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x175, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(49, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 233) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 233, x);
		e[49] = ec - (1 ^ x);
		counter[49] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(72, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 373) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 373, x);
		e[72] = ec - (1 ^ x);
		counter[72] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat33(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[24] = {73, 68, 69, 53, 54, 55, 34, 33, 11, 64, 27, 29, 31, 12, 32, 7, 6, 8, 4, 5, 3, 2, 0, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0xbf71f72047c78c69, 0xf461d6c06e45606, 0x38f6c8768fc06e94, 0x410e1e7590b9a11c, 0xdae3981d80dbcdc1, 0x283805c00266, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 24; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[8];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xc06ceccb49333d87, 0x41, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x604d80195, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xed47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x179, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[7] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 8, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[7];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x8b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x42f587, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x86f9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x13d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(64, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 317) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 317, x);
		e[64] = ec - (1 ^ x);
		counter[64] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x47210c23819, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x9eab, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(55, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 269) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 269, x);
		e[55] = ec - (1 ^ x);
		counter[55] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x453c655df, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(54, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 263) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 263, x);
		e[54] = ec - (1 ^ x);
		counter[54] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x44f76df, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(53, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 257) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 257, x);
		e[53] = ec - (1 ^ x);
		counter[53] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x3203f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(69, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 353) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 353, x);
		e[69] = ec - (1 ^ x);
		counter[69] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x24b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(68, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 349) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 349, x);
		e[68] = ec - (1 ^ x);
		counter[68] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(73, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 587) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 587, x);
		e[73] = ec - (1 ^ x);
		counter[73] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat41(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[27] = {71, 47, 46, 45, 52, 48, 62, 63, 65, 23, 67, 21, 20, 8, 64, 16, 22, 4, 35, 19, 12, 17, 68, 24, 41, 33, 34};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x5938431cccda460b, 0x5536ec1a55579cf0, 0x2fff66eeaff14754, 0x1fceb3436b97231a, 0x1eb23eb5af1d20b0, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 27; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[3];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x629f6989beb01b31, 0x850d0bc8b56ce28b, 0x603329c, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x3bcadfbb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x66bb0f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x89b1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x15d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(68, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 349) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 349, x);
		e[68] = ec - (1 ^ x);
		counter[68] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x74bc763e86672951, 0x24a7860ebe9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xd15d3599, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xc43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x155622d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x4b89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x59, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x13d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(64, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 317) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 317, x);
		e[64] = ec - (1 ^ x);
		counter[64] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x7be8fd1de28af02b, 0x951e, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x7081, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x15b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(67, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 347) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 347, x);
		e[67] = ec - (1 ^ x);
		counter[67] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x3060cbadc9af4841, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x14b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(65, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 331) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 331, x);
		e[65] = ec - (1 ^ x);
		counter[65] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xf8f2bf1fb96749, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(63, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 313) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 313, x);
		e[63] = ec - (1 ^ x);
		counter[63] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xccec09281b7f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(62, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 311) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 311, x);
		e[62] = ec - (1 ^ x);
		counter[62] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xe51547b893, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(48, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 229) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 229, x);
		e[48] = ec - (1 ^ x);
		counter[48] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xe9a58349, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(52, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 251) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 251, x);
		e[52] = ec - (1 ^ x);
		counter[52] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11b79f3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(45, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 211) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 211, x);
		e[45] = ec - (1 ^ x);
		counter[45] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1456d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(46, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 223) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 223, x);
		e[46] = ec - (1 ^ x);
		counter[46] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x16f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(47, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 227) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 227, x);
		e[47] = ec - (1 ^ x);
		counter[47] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(71, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 367) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 367, x);
		e[71] = ec - (1 ^ x);
		counter[71] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat42(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[21] = {73, 49, 58, 59, 56, 26, 53, 25, 36, 69, 15, 29, 39, 40, 70, 28, 43, 44, 37, 11, 18};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x30cdf993fb507231, 0x4d80f6e63c84b637, 0xbd2e0efba145e14d, 0xeab828e6ebd38ffc, 0x694cf5837a196436, 0x19c2a8079, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 21; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[3];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x6a881af611073d5b, 0x5bfe8184c55, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x5ecac945, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x1abf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0xa7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x79f193, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x9e77, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x167, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(70, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 359) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 359, x);
		e[70] = ec - (1 ^ x);
		counter[70] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc50bb6ff4143ec11, 0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1c386ddf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x285c25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x515b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x161, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(69, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 353) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 353, x);
		e[69] = ec - (1 ^ x);
		counter[69] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11ed7177acc8d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x6767, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x101, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(53, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 257) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 257, x);
		e[53] = ec - (1 ^ x);
		counter[53] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x288494ea9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x10f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(56, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 271) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 271, x);
		e[56] = ec - (1 ^ x);
		counter[56] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x24a6f8b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(59, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 283) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 283, x);
		e[59] = ec - (1 ^ x);
		counter[59] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x21643, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(58, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 281) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 281, x);
		e[58] = ec - (1 ^ x);
		counter[58] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x24b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(49, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 233) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 233, x);
		e[49] = ec - (1 ^ x);
		counter[49] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(73, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 587) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 587, x);
		e[73] = ec - (1 ^ x);
		counter[73] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat43(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[26] = {72, 50, 51, 57, 60, 54, 61, 66, 42, 38, 13, 55, 27, 32, 31, 14, 30, 9, 6, 10, 7, 5, 3, 2, 0, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x47d0dc5ff55c3753, 0x5db24e7101bcf938, 0x74b66000124ae2d0, 0xf248aa7eac37c85b, 0xef345a70f9af3d93, 0xd14cd0b, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 26; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[8];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xb965882316a7d2bd, 0x3d968d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x6e3e267f7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x12d67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x353, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[7] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 8, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[7];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3e3063, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x7289, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(55, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 269) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 269, x);
		e[55] = ec - (1 ^ x);
		counter[55] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1f42dd62ec6b08f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xfe11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xad, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x151, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(66, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 337) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 337, x);
		e[66] = ec - (1 ^ x);
		counter[66] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1a11660e1cb35, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(61, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 307) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 307, x);
		e[61] = ec - (1 ^ x);
		counter[61] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x195fc798ee3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(54, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 263) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 263, x);
		e[54] = ec - (1 ^ x);
		counter[54] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x162b7e567, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(60, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 293) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 293, x);
		e[60] = ec - (1 ^ x);
		counter[60] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x147d38b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(57, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 277) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 277, x);
		e[57] = ec - (1 ^ x);
		counter[57] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x15c3b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(51, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 241) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 241, x);
		e[51] = ec - (1 ^ x);
		counter[51] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x175, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(50, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 239) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 239, x);
		e[50] = ec - (1 ^ x);
		counter[50] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(72, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 373) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 373, x);
		e[72] = ec - (1 ^ x);
		counter[72] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat51(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[27] = {72, 46, 45, 47, 48, 50, 56, 69, 68, 28, 66, 15, 24, 31, 63, 25, 20, 9, 32, 67, 16, 21, 22, 37, 13, 33, 17};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x2fd271478e5ed8d5, 0xd88a116a49058417, 0x7ba513dfc2dacbcf, 0x8b4256a5eebcf183, 0x53eaa3be73ec98f, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 27; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[3];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x170b91cb3c8e80e1, 0x28a87024e3320553, 0x1c83, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x11dbf43e47d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x951dfb5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xa7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1acebd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x52af, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x15b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(67, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 347) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 347, x);
		e[67] = ec - (1 ^ x);
		counter[67] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x56dbd1bd175e941d, 0xb2896648070, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4b4bb5f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x7def, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x139, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(63, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 313) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 313, x);
		e[63] = ec - (1 ^ x);
		counter[63] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xbac68785667c7bf3, 0xae32, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1ea477, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x4dab, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x151, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(66, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 337) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 337, x);
		e[66] = ec - (1 ^ x);
		counter[66] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x217b0fd5fc19edff, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x15d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(68, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 349) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 349, x);
		e[68] = ec - (1 ^ x);
		counter[68] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xd1ef5d6ef24b5f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(69, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 353) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 353, x);
		e[69] = ec - (1 ^ x);
		counter[69] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc650a3d570b1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(56, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 271) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 271, x);
		e[56] = ec - (1 ^ x);
		counter[56] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xd46bcc685f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(50, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 239) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 239, x);
		e[50] = ec - (1 ^ x);
		counter[50] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xed7763f3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(48, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 229) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 229, x);
		e[48] = ec - (1 ^ x);
		counter[48] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x10bcdb1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(47, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 227) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 227, x);
		e[47] = ec - (1 ^ x);
		counter[47] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x144eb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(45, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 211) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 211, x);
		e[45] = ec - (1 ^ x);
		counter[45] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x175, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(46, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 223) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 223, x);
		e[46] = ec - (1 ^ x);
		counter[46] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(72, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 373) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 373, x);
		e[72] = ec - (1 ^ x);
		counter[72] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat52(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[21] = {71, 52, 57, 59, 60, 26, 64, 27, 36, 62, 23, 38, 39, 61, 19, 42, 43, 44, 18, 10, 7};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x94f923db2f6e435, 0xf66b2ed0ad58f2fa, 0x9e709a4bbf288218, 0x4c04d11ab37742a9, 0x2a5cc4a735695f55, 0xbe50e119d, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 21; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[4];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x880f508e0cff91a7, 0xf08c61a4e, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x32c9d737, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3731, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0xc7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x41ffcb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x578b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x133, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(61, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 307) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 307, x);
		e[61] = ec - (1 ^ x);
		counter[61] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1f0478e68f2d6997, 0x45, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4fa24b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x75d7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x137, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(62, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 311) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 311, x);
		e[62] = ec - (1 ^ x);
		counter[62] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xcde68c2820a5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x86f9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x13d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(64, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 317) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 317, x);
		e[64] = ec - (1 ^ x);
		counter[64] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1ae69e8c3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x125, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(60, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 293) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 293, x);
		e[60] = ec - (1 ^ x);
		counter[60] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1855979, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(59, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 283) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 283, x);
		e[59] = ec - (1 ^ x);
		counter[59] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x167d5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(57, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 277) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 277, x);
		e[57] = ec - (1 ^ x);
		counter[57] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x16f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(52, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 251) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 251, x);
		e[52] = ec - (1 ^ x);
		counter[52] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(71, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 367) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 367, x);
		e[71] = ec - (1 ^ x);
		counter[71] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat53(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[25] = {70, 49, 51, 65, 53, 54, 55, 35, 30, 14, 58, 29, 40, 41, 12, 34, 11, 6, 8, 5, 4, 3, 2, 0, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0xea0634a2361a2943, 0xe7a6a9475422032a, 0xbacaf6e56ec177d2, 0xb2eda43defc47706, 0xa8fd30a446901db6, 0x17c31087d4, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 25; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[8];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x156f9f510ce3bcf5, 0x1c01, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xc5a102c4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x1cb7d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1ed, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[7] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 8, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[7];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x628fd3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xbf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x8b67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x119, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(58, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 281) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 281, x);
		e[58] = ec - (1 ^ x);
		counter[58] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x19a307da13e83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xa4f9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(55, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 269) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 269, x);
		e[55] = ec - (1 ^ x);
		counter[55] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x18f459683a5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(54, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 263) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 263, x);
		e[54] = ec - (1 ^ x);
		counter[54] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x18db7dea5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(53, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 257) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 257, x);
		e[53] = ec - (1 ^ x);
		counter[53] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x13399cf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(65, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 331) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 331, x);
		e[65] = ec - (1 ^ x);
		counter[65] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x146bf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(51, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 241) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 241, x);
		e[51] = ec - (1 ^ x);
		counter[51] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x167, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(49, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 233) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 233, x);
		e[49] = ec - (1 ^ x);
		counter[49] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(70, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 359) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 359, x);
		e[70] = ec - (1 ^ x);
		counter[70] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat61(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[23] = {72, 46, 50, 51, 41, 61, 49, 15, 26, 58, 23, 24, 44, 14, 21, 22, 38, 45, 19, 20, 40, 35, 13};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x192ca3e788501745, 0x35d061387fdc94a5, 0x2bb4803d9b8313d6, 0xaf39c515f0082aa1, 0x1ab7f6e25572e4d6, 0x17dc026, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 23; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[3];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x52c1903391166639, 0xb89146089ca4e64, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xd20b5c9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x9d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x129145, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3c2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xd3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(45, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 211) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 211, x);
		e[45] = ec - (1 ^ x);
		counter[45] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x7c4e573a0ba2d51d, 0x3ad013e, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4a4d2a1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xd5b89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2933, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xc7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x66b8b668fb7e66e1, 0x16, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x6a79, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x119, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(58, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 281) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 281, x);
		e[58] = ec - (1 ^ x);
		counter[58] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xff8174115811, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x35b3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xe9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(49, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 233) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 233, x);
		e[49] = ec - (1 ^ x);
		counter[49] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xd50f6349ab, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(61, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 307) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 307, x);
		e[61] = ec - (1 ^ x);
		counter[61] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11d914615, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x12f5765, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(51, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 241) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 241, x);
		e[51] = ec - (1 ^ x);
		counter[51] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x144eb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(50, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 239) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 239, x);
		e[50] = ec - (1 ^ x);
		counter[50] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x175, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(46, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 223) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 223, x);
		e[46] = ec - (1 ^ x);
		counter[46] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(72, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 373) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 373, x);
		e[72] = ec - (1 ^ x);
		counter[72] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat62(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[22] = {68, 42, 43, 52, 54, 56, 29, 57, 25, 32, 59, 17, 36, 37, 60, 16, 39, 18, 10, 30, 7, 6};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x4cbb5ff0638750f9, 0x4452648b5ce99b4e, 0xec2b5c9485fe2258, 0x94a93d01189bc1db, 0xc1e4b349d8ce85de, 0x4e39da907, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 22; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[4];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xf9a2f1bb49978fa1, 0x30ae08d7a90, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1f4f21a29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x30d123, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x45d1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x125, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(60, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 293) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 293, x);
		e[60] = ec - (1 ^ x);
		counter[60] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x521d9949fe3b5d4d, 0x1951, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x2f28d3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x4a11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x11b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(59, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 283) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 283, x);
		e[59] = ec - (1 ^ x);
		counter[59] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x6b1b001808519d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x6f73, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x115, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(57, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 277) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 277, x);
		e[57] = ec - (1 ^ x);
		counter[57] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xcbf2950bed, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x10f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(56, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 271) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 271, x);
		e[56] = ec - (1 ^ x);
		counter[56] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc684f26b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(54, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 263) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 263, x);
		e[54] = ec - (1 ^ x);
		counter[54] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xca7951, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(52, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 251) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 251, x);
		e[52] = ec - (1 ^ x);
		counter[52] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1071d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x15d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(68, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 349) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 349, x);
		e[68] = ec - (1 ^ x);
		counter[68] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat63(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[22] = {64, 47, 48, 53, 62, 63, 31, 27, 9, 55, 28, 34, 12, 33, 11, 8, 4, 5, 3, 2, 0, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x9ae0cff5494cc87, 0xfa465edf8da83b04, 0x9595966b7deb7ec5, 0xfa766d3d7060010a, 0xa816325d9c7eac0d, 0x51466a68342029a, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 22; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[8];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x586c72b9fbe0792f, 0xa, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xbc397a9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x17dd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x179, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[7] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 8, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[7];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x76bd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(55, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 269) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 269, x);
		e[55] = ec - (1 ^ x);
		counter[55] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x132a8708e95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xa781, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x6d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x139, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(63, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 313) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 313, x);
		e[63] = ec - (1 ^ x);
		counter[63] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xfc6d0493, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(62, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 311) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 311, x);
		e[62] = ec - (1 ^ x);
		counter[62] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xfb7193, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(53, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 257) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 257, x);
		e[53] = ec - (1 ^ x);
		counter[53] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11917, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(48, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 229) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 229, x);
		e[48] = ec - (1 ^ x);
		counter[48] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x13d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(47, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 227) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 227, x);
		e[47] = ec - (1 ^ x);
		counter[47] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(64, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 317) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 317, x);
		e[64] = ec - (1 ^ x);
		counter[64] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat71(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[39] = {52, 42, 48, 37, 38, 39, 25, 51, 26, 20, 49, 15, 16, 2, 19, 29, 8, 9, 6, 40, 21, 22, 23, 11, 24, 17, 47, 27, 28, 30, 12, 7, 3, 33, 14, 13, 4, 5, 10};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x99b09d98c0ec0fe1, 0x156d4d0f9fb502dc, 0xf0e0e117c1b965f2, 0xded7b270b3f05069, 0x939, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 39; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[5];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x776565109da188cb, 0x89b156dff938d77d, 0x113801371fc, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x39fc0228bdb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x1ed9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2893, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x2f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2aa9b7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x60a7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xe3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(47, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 227) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 227, x);
		e[47] = ec - (1 ^ x);
		counter[47] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x55f7fb89b2bbb54b, 0x86b1482fcefe3c2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x13cf1a20f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x65, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1466d7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3aaf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xb5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x92159a820faa6da5, 0x28c, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x198a8259, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xe63, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x7f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x599191, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x35b3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xe9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(49, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 233) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 233, x);
		e[49] = ec - (1 ^ x);
		counter[49] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x14fe453b05ee31, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x64bb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xf1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(51, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 241) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 241, x);
		e[51] = ec - (1 ^ x);
		counter[51] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x4a9f55c13d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xb3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x6e6c84d1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xa945c7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xbd3b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(48, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 229) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 229, x);
		e[48] = ec - (1 ^ x);
		counter[48] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xfb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(52, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 251) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 251, x);
		e[52] = ec - (1 ^ x);
		counter[52] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat72(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[17] = {57, 41, 0, 1, 43, 44, 45, 46, 56, 32, 54, 35, 53, 36, 31, 34, 18};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x6e80e2c52d6e1489, 0x13954dde593a80be, 0x9b177a8227248fe8, 0xf2a90f95b93f7a2e, 0xe54d2323c17a3663, 0x3e62efdcaa2d4ad7, 0xf1, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 17; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[3];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xdbb7339ad043eb65, 0x1e228, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x57923b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xa3a3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x101, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(53, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 257) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 257, x);
		e[53] = ec - (1 ^ x);
		counter[53] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xfd456b22af9c680f, 0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x107, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(54, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 263) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 263, x);
		e[54] = ec - (1 ^ x);
		counter[54] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x53368783a0223, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x10f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(56, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 271) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 271, x);
		e[56] = ec - (1 ^ x);
		counter[56] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x5f86ebed33d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(46, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 223) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 223, x);
		e[46] = ec - (1 ^ x);
		counter[46] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x73e66d1af, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(45, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 211) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 211, x);
		e[45] = ec - (1 ^ x);
		counter[45] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x9518fd9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc1c05, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x115, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xbf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(57, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 277) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 277, x);
		e[57] = ec - (1 ^ x);
		counter[57] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat81(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[45] = {44, 29, 40, 42, 43, 41, 15, 39, 20, 26, 35, 16, 27, 28, 8, 36, 24, 19, 31, 37, 21, 22, 23, 17, 25, 9, 13, 12, 4, 3, 38, 32, 33, 34, 18, 30, 14, 6, 11, 7, 10, 5, 2, 0, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x361d8dba571a5213, 0x1820585358c06619, 0x5e5c0bf97a2a0127, 0x31797e52c32f, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 45; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[8];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x9161cc8845f84c59, 0x337b41346488c650, 0xf7f8c49d06c4cf8d, 0x1, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x8f19d96bb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x2034d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x3af, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[7] = R;
	{u512 coef = { .c = {0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 8, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[7];
	ec = lookup(0, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 3) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 3, x);
		e[0] = ec - (1 ^ x);
		counter[0] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x35, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x36ac1b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x5def, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xad, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x827d956b05e0aa75, 0xe896fec5882d51f6, 0x4, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1dddf2077, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x24a37, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0xc79, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x12d2dd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3625, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xa7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xfdc5e5e2edf2c24b, 0x800b664bfa, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x125687, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x404f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xa3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc38d438d28beeffb, 0xa0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xfedb5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x71, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2569, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x9d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x6f68fad9fc085, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x373d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xb3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x287ea9ea81, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xbf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(41, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 191) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 191, x);
		e[41] = ec - (1 ^ x);
		counter[41] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x349f668d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(43, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 197) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 197, x);
		e[43] = ec - (1 ^ x);
		counter[43] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x45cccd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(42, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 193) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 193, x);
		e[42] = ec - (1 ^ x);
		counter[42] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x62b9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xc7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(44, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 199) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 199, x);
		e[44] = ec - (1 ^ x);
		counter[44] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat91(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[40] = {40, 29, 37, 35, 15, 39, 16, 2, 27, 8, 38, 24, 25, 9, 34, 20, 28, 19, 13, 33, 21, 22, 23, 14, 26, 12, 18, 36, 30, 31, 17, 32, 11, 5, 10, 6, 7, 4, 3, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0xb24a588f2e1e7535, 0xc2b8a8bbd7be633, 0x71eac4fe166e8ec0, 0xc965acdf4a302dd3, 0x31f0, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 40; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[7];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xb97665551f5e7179, 0x930c661c30485d6a, 0xaf3befd72b7, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xbaeb5d3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x17a73, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2bf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x8b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x5369, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xa3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xe7707d9703bd92bd, 0x1aaf067d0d95a94, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x15145fcc3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x11f9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x6b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x10cb77, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x304f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x7d9bae5a87372aa3, 0x18c7b2e, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x149189, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2e99, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x773a36ed9731363b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4441, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xad, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x8477c6503, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x12a91, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x6d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xb3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(39, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 179) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 179, x);
		e[39] = ec - (1 ^ x);
		counter[39] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x3a936d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x9d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x59cb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xb5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(40, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 181) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 181, x);
		e[40] = ec - (1 ^ x);
		counter[40] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat101(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[38] = {38, 29, 32, 37, 31, 22, 35, 15, 23, 33, 27, 26, 30, 34, 19, 20, 18, 10, 14, 8, 3, 36, 16, 24, 25, 11, 28, 13, 12, 4, 2, 21, 9, 7, 17, 6, 5, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x2f1df0adc0b38e9b, 0x3a832b51d38308a1, 0x2f1918f8f2044344, 0x8ae1773e2777ab2e, 0x18b068a0, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 38; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[6];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x91f855cebde623e9, 0x9385c97b9add772e, 0x90, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1391de1094f196b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xe72b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x13, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x53, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xfcc746ad, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x14bf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x71, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xf52d3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x26d7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xa3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(36, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 163) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 163, x);
		e[36] = ec - (1 ^ x);
		counter[36] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xadfda0b972db6369, 0x41c29de, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xd49a1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xa43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x35, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2b0f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(34, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 151) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 151, x);
		e[34] = ec - (1 ^ x);
		counter[34] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x4d8aacc2df6e46b9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1a843b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3f71, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x95, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x5a7e0a164f7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x242f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x9d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1e662d97, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x2e9991, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(37, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 167) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 167, x);
		e[37] = ec - (1 ^ x);
		counter[37] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x55d3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xad, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(38, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 173) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 173, x);
		e[38] = ec - (1 ^ x);
		counter[38] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat111(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[34] = {35, 28, 31, 33, 27, 22, 16, 9, 32, 15, 25, 20, 29, 14, 21, 18, 8, 4, 19, 13, 11, 6, 5, 30, 26, 23, 17, 24, 12, 10, 7, 3, 2, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x5c049208954d3c1d, 0xd6034aa8aa77786, 0x2349400abaa6233, 0x340d7a8dbdb1d741, 0x4166f76fcb004d2, 0x0, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 34; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[7];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x3dbd425743fd4af5, 0x373aa31bc1afba04, 0x375, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x56e05a3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x10f7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[6] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 7, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[6];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x65, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x36c1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x83, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xccc58b1c76bc6ca7, 0x5de0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x37b53ef3f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xd67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x13, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x88651, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1a4b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x7f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x179a2c39bca737, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xce39f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x2009, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x8b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(32, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 139) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 139, x);
		e[32] = ec - (1 ^ x);
		counter[32] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1595ecd9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x25e5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x6d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x251635, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(33, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 149) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 149, x);
		e[33] = ec - (1 ^ x);
		counter[33] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x454d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(31, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 137) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 137, x);
		e[31] = ec - (1 ^ x);
		counter[31] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x9d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(35, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 157) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 157, x);
		e[35] = ec - (1 ^ x);
		counter[35] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat121(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[30] = {29, 27, 30, 24, 15, 4, 28, 16, 2, 18, 26, 13, 14, 21, 9, 17, 8, 7, 3, 1, 25, 22, 23, 20, 11, 19, 10, 12, 6, 5};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x1afe06066b09a2d7, 0x3ffe47a43addb19c, 0x6e35d789027d887a, 0x9630687aebc11515, 0xfb7ab61550f154f4, 0x6c8937, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 30; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[5];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x1e31f87357530955, 0x27ff26b09273, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xabab1369, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xa8d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x13, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xd916f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x23cf, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x6adca99044e0d57, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x28e07915, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x797, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x41129, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x53, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x13a5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x6b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x20b5553703, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xbc7b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x71, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1babb9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x65, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x3613, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(30, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 131) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 131, x);
		e[30] = ec - (1 ^ x);
		counter[30] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x7f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(29, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 127) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 127, x);
		e[29] = ec - (1 ^ x);
		counter[29] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat131(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[28] = {28, 23, 26, 27, 16, 4, 25, 14, 12, 8, 24, 13, 18, 15, 9, 7, 3, 2, 1, 22, 21, 20, 19, 17, 10, 11, 6, 5};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x2e8d8f1a313bae7b, 0xd035dc05997f8e06, 0x61184306d82e70d0, 0x87f5fc6b431731bb, 0x35cb474d53be1b44, 0x1b8d919b32, 0x0, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 28; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[6];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xde9e8480374cd3ad, 0x2f8bc2dabb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x28a097d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x9af, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x13, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x8e795, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x1cdb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x59, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x993bbf8bfc5313b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x5248d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x725, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x3b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x128b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x65, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(24, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 101) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 101, x);
		e[24] = ec - (1 ^ x);
		counter[24] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x179a759ea7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x1553, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x67, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(25, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 103) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 103, x);
		e[25] = ec - (1 ^ x);
		counter[25] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11e55b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x6d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(27, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 109) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 109, x);
		e[27] = ec - (1 ^ x);
		counter[27] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x2ad1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(26, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 107) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 107, x);
		e[26] = ec - (1 ^ x);
		counter[26] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x71, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(28, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 113) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 113, x);
		e[28] = ec - (1 ^ x);
		counter[28] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat141(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[23] = {23, 15, 22, 16, 9, 20, 8, 14, 21, 11, 12, 6, 13, 7, 4, 2, 1, 19, 18, 17, 10, 5, 3};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x868446907f4395cf, 0xbee3acd3f5e954a4, 0xb82e0674b1b9a334, 0x85ff827d8db1e96d, 0x16a265b8b2fcdba6, 0xf443f9b2b4c64a76, 0x57, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 23; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[6];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xc0971771beaa8a53, 0x844a4a, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x143f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x6a5debbcb9c7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x2a6c5b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x2f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xd4b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x53, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x3968e309, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x8f3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x7c5a3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x165b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat151(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[23] = {23, 16, 19, 20, 12, 22, 13, 21, 9, 7, 11, 4, 3, 18, 15, 17, 6, 14, 8, 2, 1, 10, 5};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x868446907f4395cf, 0xbee3acd3f5e954a4, 0xb82e0674b1b9a334, 0x85ff827d8db1e96d, 0x16a265b8b2fcdba6, 0xf443f9b2b4c64a76, 0x57, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 23; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[5];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x71f458d49bb56c33, 0x73, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x515e75, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0xd223, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x35, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x105d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x5950d448dff, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xe72b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x29, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x53, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x1f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x57754699, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x59, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x69745, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x171d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x61, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(23, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 97) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 97, x);
		e[23] = ec - (1 ^ x);
		counter[23] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat161(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[22] = {22, 20, 21, 18, 11, 19, 14, 17, 16, 9, 12, 13, 8, 4, 2, 15, 10, 7, 6, 5, 3, 1};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0xf81ebcc0389bc36f, 0x54447c502d691256, 0xc97072375756d6fc, 0xc5d07192b0697292, 0x93888afbd1cf3a18, 0x8dc19cb67f2236be, 0x2153, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 22; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[6];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x19082069ef93d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x259247221, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x887, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1b5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0xb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x3b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x13da1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x2f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x763, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x6541574b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xf1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x8e795, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1b77, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x59, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(22, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 89) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 89, x);
		e[22] = ec - (1 ^ x);
		counter[22] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat171(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[21] = {21, 15, 20, 13, 19, 8, 7, 16, 9, 11, 12, 4, 2, 1, 18, 10, 17, 14, 6, 5, 3};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x42af9ed3ae26f197, 0x4bcf37dfc987603c, 0x817b53d5d30bdb9, 0xc5777bff54a8d508, 0x4a78518bf10b329c, 0x484f7b7232e50841, 0xb960c, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 21; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[5];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0xcd23676830e178ad, 0x48, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x2af89, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x35, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x13, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0xa43, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x47, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0xce219b2113, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x12edb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(1, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 5) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 5, x);
		e[1] = ec - (1 ^ x);
		counter[1] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x763, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x11571a1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x1d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x1321, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x53, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(21, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 83) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 83, x);
		e[21] = ec - (1 ^ x);
		counter[21] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void strat181(proj *A, int8_t *e, int8_t *counter, unsigned int *isog_counter){

	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};
	int8_t ec = 0;
	uint8_t bc, x;
	proj P, S, R;
	u512 k;

	// Prime indices used in strategy:
	int actv_ind[19] = {20, 18, 19, 14, 15, 13, 17, 16, 9, 11, 12, 6, 5, 10, 8, 7, 4, 2, 3};

	// Choose a random point P with Elligator. When on the base curve, we choose a point of full order
	if(memcmp(&A->x, &fp_0, sizeof(fp))) {
		elligator(&P.x, &A->x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;
	} else {
		fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
		P.z = fp_1;
	}

	// Multiply P by any primes not used in strategy:
	{u512 coef = { .c = {0x1ab279275121a3c9, 0xe4eb93c7b27501b0, 0x1e6ece7a12038f61, 0x1cb202ea3db15805, 0xb90c31dbc1270c24, 0x38d91e20814861d7, 0x12c83de9, 0x0}};
	xMUL(&P, A, &P, &coef);}

	// Multiply P by any primes for which a real isogeny wont be constructed
	S = P;
	for (uint8_t i = 0; i < 19; i++){
		ec = lookup(actv_ind[i], e);  //check in constant-time if normal or dummy isogeny must be computed
		bc = isequal(ec, 0);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
		u512_set(&k, primes[actv_ind[i]]);
		xMUL(&S, A, &S, &k);
		fp_cswap(&P.x, &S.x, bc);
		fp_cswap(&P.z, &S.z, bc);
	}
	xDBL(&P, A, &P);		// Remove 4-torsion from P
	xDBL(&P, A, &P);


	// Initialize array to store points for isogeny evaluation
	proj pts[6];
	pts[0] = P;
	R = pts[0];
	{u512 coef = { .c = {0x3aabf76b25b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x402f179b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x431, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x17, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[4] = R;
	{u512 coef = { .c = {0xd, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[5] = R;
	{u512 coef = { .c = {0x7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(3, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 11) == 0 ));
	if (x < 2){
		xISOG(A, pts, 6, &R, 11, x);
		e[3] = ec - (1 ^ x);
		counter[3] -= 1;
		(*isog_counter)++;
	}
	R = pts[5];
	ec = lookup(2, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 7) == 0 ));
	if (x < 2){
		xISOG(A, pts, 5, &R, 7, x);
		e[2] = ec - (1 ^ x);
		counter[2] -= 1;
		(*isog_counter)++;
	}
	R = pts[4];
	ec = lookup(4, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 13) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 13, x);
		e[4] = ec - (1 ^ x);
		counter[4] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(7, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 23) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 23, x);
		e[7] = ec - (1 ^ x);
		counter[7] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	{u512 coef = { .c = {0x25, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(8, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 29) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 29, x);
		e[8] = ec - (1 ^ x);
		counter[8] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(10, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 37) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 37, x);
		e[10] = ec - (1 ^ x);
		counter[10] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x12edb, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[2] = R;
	{u512 coef = { .c = {0x2b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[3] = R;
	{u512 coef = { .c = {0x13, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(5, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 17) == 0 ));
	if (x < 2){
		xISOG(A, pts, 4, &R, 17, x);
		e[5] = ec - (1 ^ x);
		counter[5] -= 1;
		(*isog_counter)++;
	}
	R = pts[3];
	ec = lookup(6, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 19) == 0 ));
	if (x < 2){
		xISOG(A, pts, 3, &R, 19, x);
		e[6] = ec - (1 ^ x);
		counter[6] -= 1;
		(*isog_counter)++;
	}
	R = pts[2];
	ec = lookup(12, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 43) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 43, x);
		e[12] = ec - (1 ^ x);
		counter[12] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x763, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(11, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 41) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 41, x);
		e[11] = ec - (1 ^ x);
		counter[11] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3d, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(9, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 31) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 31, x);
		e[9] = ec - (1 ^ x);
		counter[9] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(16, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 61) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 61, x);
		e[16] = ec - (1 ^ x);
		counter[16] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x14b2265, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0xad5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(17, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 67) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 67, x);
		e[17] = ec - (1 ^ x);
		counter[17] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	{u512 coef = { .c = {0x3b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(13, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 47) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 47, x);
		e[13] = ec - (1 ^ x);
		counter[13] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(15, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 59) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 59, x);
		e[15] = ec - (1 ^ x);
		counter[15] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x15e9, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	pts[1] = R;
	{u512 coef = { .c = {0x49, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(14, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 53) == 0 ));
	if (x < 2){
		xISOG(A, pts, 2, &R, 53, x);
		e[14] = ec - (1 ^ x);
		counter[14] -= 1;
		(*isog_counter)++;
	}
	R = pts[1];
	ec = lookup(19, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 73) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 73, x);
		e[19] = ec - (1 ^ x);
		counter[19] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	{u512 coef = { .c = {0x4f, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
	xMUL(&R, A, &R, &coef);}
	ec = lookup(18, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 71) == 0 ));
	if (x < 2){
		xISOG(A, pts, 1, &R, 71, x);
		e[18] = ec - (1 ^ x);
		counter[18] -= 1;
		(*isog_counter)++;
	}
	R = pts[0];
	ec = lookup(20, e);
	bc = 1 ^ isequal(ec, 0);
	x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * 79) == 0 ));
	if (x < 2){
		lastxISOG(A, &R, 79, x);
		e[20] = ec - (1 ^ x);
		counter[20] -= 1;
		(*isog_counter)++;
	}

	// Affinize the curve coefficients
	fp_inv(&A->z);
	fp_mul2(&A->x, &A->z);
	A->z = fp_1;

}


void action(public_key *out, public_key const *in, private_key const *priv,int8_t const *max_exponent) {
	int8_t e[num_primes] = {0};
	unsigned int isog_counter = 0;
	int8_t counter[num_primes] = {0};
	uint8_t x;

	memcpy(e, priv->e, sizeof(priv->e));			// Copy of the private key
	memcpy(counter, max_exponent, sizeof(counter));	// Counts how many isogenies have been constructed for each prime

	proj A = { in->A, fp_1 };	 // Projectivize the curve coefficient

	strat11(&A, e, counter, &isog_counter);
	strat12(&A, e, counter, &isog_counter);
	strat13(&A, e, counter, &isog_counter);
	strat21(&A, e, counter, &isog_counter);
	strat22(&A, e, counter, &isog_counter);
	strat23(&A, e, counter, &isog_counter);
	strat31(&A, e, counter, &isog_counter);
	strat32(&A, e, counter, &isog_counter);
	strat33(&A, e, counter, &isog_counter);
	strat41(&A, e, counter, &isog_counter);
	strat42(&A, e, counter, &isog_counter);
	strat43(&A, e, counter, &isog_counter);
	strat51(&A, e, counter, &isog_counter);
	strat52(&A, e, counter, &isog_counter);
	strat53(&A, e, counter, &isog_counter);
	strat61(&A, e, counter, &isog_counter);
	strat62(&A, e, counter, &isog_counter);
	strat63(&A, e, counter, &isog_counter);
	strat71(&A, e, counter, &isog_counter);
	strat72(&A, e, counter, &isog_counter);
	strat81(&A, e, counter, &isog_counter);
	strat91(&A, e, counter, &isog_counter);
	strat101(&A, e, counter, &isog_counter);
	strat111(&A, e, counter, &isog_counter);
	strat121(&A, e, counter, &isog_counter);
	strat131(&A, e, counter, &isog_counter);
	strat141(&A, e, counter, &isog_counter);
	strat151(&A, e, counter, &isog_counter);
	strat161(&A, e, counter, &isog_counter);
	strat171(&A, e, counter, &isog_counter);
	strat181(&A, e, counter, &isog_counter);


	int8_t ec = 0;
	uint8_t bc;
	proj P, S, R;
	u512 k;
	// Now finish any isogeny constructions which failed above
	while (isog_counter != 815){

		// Choose a random point P with Elligator.
		elligator(&P.x, &A.x, true, elligator_index);
		elligator_index = (elligator_index + 1)%9;
		P.z = fp_1;

		// Multiply P by any primes for which a real isogeny wont be constructed
		S = P;
		for (uint8_t i = 0; i < 74; i++){
			ec = lookup(i, e);  //check in constant-time if normal or dummy isogeny must be computed
			bc = isequal(ec, 0);
			fp_cswap(&P.x, &S.x, !bc);
			fp_cswap(&P.z, &S.z, !bc);
			u512_set(&k, primes[i]);
			xMUL(&P, &A, &P, &k);
			fp_cswap(&P.x, &S.x, !bc);
			fp_cswap(&P.z, &S.z, !bc);
		}
		xDBL(&P, &A, &P);		// Remove 4-torsion
		xDBL(&P, &A, &P);
		for (int i = 0; i < num_primes; i++){
			if (counter[i] != 0){
				R = P;
				for (int j = i+1; j < num_primes; j++ ){
					if (counter[j] != 0){
						u512_set(&k, primes[j]);
						xMUL(&R, &A, &R, &k);
					}
				}
				ec = lookup(i, e);
				bc = 1 ^ isequal(ec, 0);
				x = 2 - 2*(bc && memcmp(&R.z, &fp_0, sizeof(fp))) - (1 ^ (bc || (int)((float)rand()/RAND_MAX * primes[i]) == 0 ));
				if (x < 2){
					xISOG(&A, &P, 1, &R, primes[i], x);
					e[i] = ec - (1 ^ x);
					counter[i] -= 1;
					isog_counter++;
				}
			}
		}
		fp_inv(&A.z);
		fp_mul2(&A.x, &A.z);
		A.z = fp_1;


	}

	out->A = A.x;

}



/* includes public-key validation. */
bool csidh(public_key *out, public_key const *in, private_key const *priv, int8_t const *max_exponent) {
	if (!validate(in)) {
		fp_random(&out->A);
		return false;
	}
	action(out, in, priv, max_exponent);

	return true;
}
