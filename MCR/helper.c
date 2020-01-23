
#include <string.h>
#include <assert.h>

#include "helper.h"
#include "rng.h"

const u512 four_sqrt_p = { { 0x85e2579c786882cf, 0x4e3433657e18da95,
		0x850ae5507965a0b3, 0xa15bc4e676475964, } };

const public_key base = { 0 }; /* A = 0 */

/* get priv[pos] in constant time  */
uint32_t lookup(size_t pos, int8_t const *priv)
{
	int b;
	int8_t r = priv[0];
	for(size_t i=1;i<num_primes;i++)
	{
		b = isequal(i, pos);
		//ISEQUAL(i, pos, b);
		//b = (uint8_t)(1-((-(i ^ pos)) >> 31));
		cmov(&r, &priv[i], b);
		//CMOV(&r, &priv[i], b);
	}
	return r;
}

/* check if a and b are equal in constant time  */
uint32_t isequal(uint32_t a, uint32_t b)
{
	//size_t i;
	uint32_t r = 0;
	unsigned char *ta = (unsigned char *)&a;
	unsigned char *tb = (unsigned char *)&b;
	r = (ta[0] ^ tb[0]) | (ta[1] ^ tb[1]) | (ta[2] ^ tb[2]) |  (ta[3] ^ tb[3]);
	r = (-r);
	r = r >> 31;
	return (int)(1-r);
}


/* decision bit b has to be either 0 or 1 */
void cmov(int8_t *r, const int8_t *a, uint32_t b)
{
	uint32_t t;
	b = -b; /* Now b is either 0 or 0xffffffff */
	t = (*r ^ *a) & b;
	*r ^= t;
}


void csidh_private(private_key *priv, const int8_t *max_exponent) {
	memset(&priv->e, 0, sizeof(priv->e));
	for (size_t i = 0; i < num_primes;) {
		int8_t buf[64];
		randombytes(buf, sizeof(buf));
		for (size_t j = 0; j < sizeof(buf); ++j) {
			if (buf[j] <= max_exponent[i] && buf[j] >= 0) {
				priv->e[i] = lookup(j, buf);
				if (++i >= num_primes)
					break;
			}
		}
	}
}

/* compute [(p+1)/l] P for all l in our list of primes. */
/* divide and conquer is much faster than doing it naively,
 * but uses more memory. */
static void cofactor_multiples(proj *P, const proj *A, size_t lower,
		size_t upper) {
	assert(lower < upper);

	if (upper - lower == 1)
		return;

	size_t mid = lower + (upper - lower + 1) / 2;

	u512 cl = u512_1, cu = u512_1;
	for (size_t i = lower; i < mid; ++i)
		u512_mul3_64(&cu, &cu, primes[i]);
	for (size_t i = mid; i < upper; ++i)
		u512_mul3_64(&cl, &cl, primes[i]);

	xMUL(&P[mid], A, &P[lower], &cu);
	xMUL(&P[lower], A, &P[lower], &cl);

	cofactor_multiples(P, A, lower, mid);
	cofactor_multiples(P, A, mid, upper);
}

/* never accepts invalid keys. */
bool validate(public_key const *in) {
	const proj A = { in->A, fp_1 };

	do {

		proj P[num_primes];
		fp_random(&P->x);
		P->z = fp_1;

		/* maximal 2-power in p+1 */
		xDBL(P, &A, P);
		xDBL(P, &A, P);

		cofactor_multiples(P, &A, 0, num_primes);

		u512 order = u512_1;

		for (size_t i = num_primes - 1; i < num_primes; --i) {

			/* we only gain information if [(p+1)/l] P is non-zero */
			if (memcmp(&P[i].z, &fp_0, sizeof(fp))) {

				u512 tmp;
				u512_set(&tmp, primes[i]);
				xMUL(&P[i], &A, &P[i], &tmp);

				if (memcmp(&P[i].z, &fp_0, sizeof(fp)))
					/* P does not have order dividing p+1. */
					return false;

				u512_mul3_64(&order, &order, primes[i]);

				if (u512_sub3(&tmp, &four_sqrt_p, &order)) /* returns borrow */
					/* order > 4 sqrt(p), hence definitely supersingular */
					return true;
			}
		}

		/* P didn't have big enough order to prove supersingularity. */
	} while (1);
}

/* compute x^3 + Ax^2 + x */
static void montgomery_rhs(fp *rhs, fp const *A, fp const *x) {
	fp tmp;
	*rhs = *x;
	fp_sq1(rhs);
	fp_mul3(&tmp, A, x);
	fp_add2(rhs, &tmp);
	fp_add2(rhs, &fp_1);
	fp_mul2(rhs, x);
}


/* generates a curve point with suitable field of definition for y-coordinate */
void elligator(fp * x, const fp *A, bool sign, uint8_t index) {

	fp legendre_symbol;
	// v = A/(u^2 − 1)
	fp_set(x, 0);
	fp_add2(x, &invs_[index]);
	fp_mul2(x, A);
	// Compute the Legendre symbol
	montgomery_rhs(&legendre_symbol, A, x);
	// Compute x as v if e = s
	if(fp_issquare(&legendre_symbol)!=sign){
		// otherwise − v − A
		fp_add2(x, A);
		fp_sub3(x, &fp_0, x);

	}

}

