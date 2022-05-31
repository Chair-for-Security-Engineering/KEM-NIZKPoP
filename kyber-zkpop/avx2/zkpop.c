#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include "zkpop.h"
#include "cbd.h"
#include "fips202.h"
#include "fips202x4.h"
#include <immintrin.h>
#include "consts.h"
#include "randombytes.h"
#include "symmetric.h"
#include "reduce.h"
#include "rejsample.h"
#include "log.h"
#include "ntt.h"

// static void printhash(const char* str, uint8_t x[ZKPOP_SYMBYTES])
// {
//   printf("%s", str);
//   for (size_t i = 0; i < ZKPOP_SYMBYTES; i++)
//   {
//     printf("%02x", x[i]);
//   }
//   printf("\n");
// }

#if ZKPOP_N >= 8192
#define BIGPARAM
#endif

#if KYBER_K == 2
#define SHAKE                 shake128
#define SHAKEx4               shake128x4
#define SHAKEx4_absorb_once   shake128x4_absorb_once
#define SHAKEx4_squeezeblocks shake128x4_squeezeblocks
#define SHAKE_init            shake128_init
#define SHAKE_absorb          shake128_absorb
#define SHAKE_absorb_once     shake128_absorb_once
#define SHAKE_finalize        shake128_finalize
#define SHAKE_squeeze         shake128_squeeze
#define SHAKE_squeezeblocks   shake128_squeezeblocks
#define SHAKE_RATE            SHAKE128_RATE
#elif KYBER_K == 3 || KYBER_K == 4
#define SHAKE                 shake256
#define SHAKEx4               shake256x4
#define SHAKEx4_absorb_once   shake256x4_absorb_once
#define SHAKEx4_squeezeblocks shake256x4_squeezeblocks
#define SHAKE_init            shake256_init
#define SHAKE_absorb          shake256_absorb
#define SHAKE_absorb_once     shake256_absorb_once
#define SHAKE_finalize        shake256_finalize
#define SHAKE_squeeze         shake256_squeeze
#define SHAKE_squeezeblocks   shake256_squeezeblocks
#define SHAKE_RATE            SHAKE256_RATE
#endif

/*
 * Macros
 */
// load byte array from zkpop buffer
#define LOADZKPOP(dest, size) do {if (zkpop_size < size) { /*printf("bad zkpop: size\n");*/ return __LINE__;} memcpy(dest, zkpop_cur, size); zkpop_cur += size; zkpop_size -= size;} while(0)

// load und unpack polyvec from zkpop buffer
#define LOADZKPOP_POLYVEC(dest) do {if (zkpop_size < KYBER_POLYVECBYTES) { /*printf("bad zkpop: size\n");*/ return __LINE__;} polyvec_frombytes(dest, zkpop_cur); zkpop_cur += KYBER_POLYVECBYTES; zkpop_size -= KYBER_POLYVECBYTES;} while(0)

// allocate array
#define CALLOC8(x, size) do {x = calloc(size, sizeof(uint8_t)); if (x == NULL) {return __LINE__;}} while(0)

// store byte array to zkpop buffer
#define STOREZKPOP(src, size) do {memcpy(zkpop_cur, src, size); zkpop_cur += size; proof_size -= size;} while(0)

// pack and store polyvec to zkpop buffer
#define STOREZKPOP_POLYVEC(src) do {polyvec_tobytes(zkpop_cur, src); zkpop_cur += KYBER_POLYVECBYTES; proof_size -= KYBER_POLYVECBYTES;} while(0)

/*
 * HELPER FUNCTIONS COPIED FROM OTHER FILES
 */

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk and the
*              public seed used to generate the matrix A.
*              The polynomial coefficients in pk are assumed to
*              lie in the invertal [0,q], i.e. pk must be reduced
*              by polyvec_reduce().
*
* Arguments:   uint8_t *r: pointer to the output serialized public key
*              polyvec *pk: pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  polyvec_tobytes(r, pk);
  memcpy(r+KYBER_POLYVECBYTES, seed, KYBER_SYMBYTES);
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk: pointer to output public-key polynomial vector
*              - uint8_t *seed: pointer to output seed to generate matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  polyvec_frombytes(pk, packedpk);
  memcpy(seed, packedpk+KYBER_POLYVECBYTES, KYBER_SYMBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key.
*              The polynomial coefficients in sk are assumed to
*              lie in the invertal [0,q], i.e. sk must be reduced
*              by polyvec_reduce().
*
* Arguments:   - uint8_t *r: pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        pack_bundle
**************************************************/
static void pack_bundle(uint8_t *r, const bundle_t *x)
{
  size_t n_x = ZKPOP_M, i=0;
  uint8_t *rptr = r;
  while (n_x >= KYBER_N)
  {
    ntttobytes_avx(rptr, &x->vec[i], qdata.vec);
    rptr += KYBER_POLYBYTES;
    n_x -= KYBER_N;
    i += KYBER_N/16;
  }
  memcpy(rptr, &x->coeffs[i*16], n_x*sizeof(int16_t));
}

/*************************************************
* Name:        unpack_bundle
**************************************************/
static void unpack_bundle(bundle_t *x, const uint8_t *r)
{
  size_t n_x = ZKPOP_M, i = 0;
  const uint8_t *rptr = r;
  while (n_x >= KYBER_N)
  {
    nttfrombytes_avx(&x->vec[i], rptr, qdata.vec);
    rptr += KYBER_POLYBYTES;
    n_x -= KYBER_N;
    i += KYBER_N/16;
  }
  memcpy(&x->coeffs[i*16], rptr, n_x*sizeof(int16_t));
}

/*************************************************
* Name:        pack_open_bundle
**************************************************/
static void pack_open_bundle(uint8_t *r, const open_bundle_t *x)
{
  size_t n_x = ZKPOP_M-ZKPOP_SIGMA, i=0;
  uint8_t *rptr = r;
  while (n_x >= KYBER_N)
  {
    ntttobytes_avx(rptr, &x->vec[i], qdata.vec);
    rptr += KYBER_POLYBYTES;
    n_x -= KYBER_N;
    i += KYBER_N/16;
  }
  memcpy(rptr, &x->coeffs[i*16], n_x*sizeof(int16_t));
}

/*************************************************
* Name:        unpack_open_bundle
**************************************************/
static void unpack_open_bundle(open_bundle_t *x, const uint8_t *r)
{
  size_t n_x = ZKPOP_M-ZKPOP_SIGMA, i = 0;
  const uint8_t *rptr = r;
  while (n_x >= KYBER_N)
  {
    nttfrombytes_avx(&x->vec[i], rptr, qdata.vec);
    rptr += KYBER_POLYBYTES;
    n_x -= KYBER_N;
    i += KYBER_N/16;
  }
  memcpy(&x->coeffs[i*16], rptr, n_x*sizeof(int16_t));
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output array
*              - unsigned int len: requested number of 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos <= buflen - 3) {  // buflen is always at least 3
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < KYBER_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#if KYBER_K == 2
static void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 0;
    buf[1].coeffs[33] = 1;
    buf[2].coeffs[32] = 1;
    buf[2].coeffs[33] = 0;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 1;
  }
  else {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 1;
    buf[1].coeffs[33] = 0;
    buf[2].coeffs[32] = 0;
    buf[2].coeffs[33] = 1;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 1;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1].coeffs);
  ctr2 = rej_uniform_avx(a[1].vec[0].coeffs, buf[2].coeffs);
  ctr3 = rej_uniform_avx(a[1].vec[1].coeffs, buf[3].coeffs);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[1].vec[0].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[1].vec[1].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[1].vec[0]);
  poly_nttunpack(&a[1].vec[1]);
}
#elif KYBER_K == 3
static void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;
  keccak_state state1x;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 0;
    buf[1].coeffs[33] = 1;
    buf[2].coeffs[32] = 0;
    buf[2].coeffs[33] = 2;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 0;
  }
  else {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 1;
    buf[1].coeffs[33] = 0;
    buf[2].coeffs[32] = 2;
    buf[2].coeffs[33] = 0;
    buf[3].coeffs[32] = 0;
    buf[3].coeffs[33] = 1;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1].coeffs);
  ctr2 = rej_uniform_avx(a[0].vec[2].coeffs, buf[2].coeffs);
  ctr3 = rej_uniform_avx(a[1].vec[0].coeffs, buf[3].coeffs);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[0].vec[2].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[1].vec[0].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[0].vec[2]);
  poly_nttunpack(&a[1].vec[0]);

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[32] = 1;
    buf[0].coeffs[33] = 1;
    buf[1].coeffs[32] = 1;
    buf[1].coeffs[33] = 2;
    buf[2].coeffs[32] = 2;
    buf[2].coeffs[33] = 0;
    buf[3].coeffs[32] = 2;
    buf[3].coeffs[33] = 1;
  }
  else {
    buf[0].coeffs[32] = 1;
    buf[0].coeffs[33] = 1;
    buf[1].coeffs[32] = 2;
    buf[1].coeffs[33] = 1;
    buf[2].coeffs[32] = 0;
    buf[2].coeffs[33] = 2;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 2;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[1].vec[1].coeffs, buf[0].coeffs);
  ctr1 = rej_uniform_avx(a[1].vec[2].coeffs, buf[1].coeffs);
  ctr2 = rej_uniform_avx(a[2].vec[0].coeffs, buf[2].coeffs);
  ctr3 = rej_uniform_avx(a[2].vec[1].coeffs, buf[3].coeffs);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[1].vec[1].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[1].vec[2].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[2].vec[0].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[2].vec[1].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[1].vec[1]);
  poly_nttunpack(&a[1].vec[2]);
  poly_nttunpack(&a[2].vec[0]);
  poly_nttunpack(&a[2].vec[1]);

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  buf[0].coeffs[32] = 2;
  buf[0].coeffs[33] = 2;
  shake128_absorb_once(&state1x, buf[0].coeffs, 34);
  shake128_squeezeblocks(buf[0].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state1x);
  ctr0 = rej_uniform_avx(a[2].vec[2].coeffs, buf[0].coeffs);
  while(ctr0 < KYBER_N) {
    shake128_squeezeblocks(buf[0].coeffs, 1, &state1x);
    ctr0 += rej_uniform(a[2].vec[2].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[2].vec[2]);
}
#elif KYBER_K == 4
static void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int i, ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  for(i=0;i<4;i++) {
    f = _mm256_loadu_si256((__m256i *)seed);
    _mm256_store_si256(buf[0].vec, f);
    _mm256_store_si256(buf[1].vec, f);
    _mm256_store_si256(buf[2].vec, f);
    _mm256_store_si256(buf[3].vec, f);

    if(transposed) {
      buf[0].coeffs[32] = i;
      buf[0].coeffs[33] = 0;
      buf[1].coeffs[32] = i;
      buf[1].coeffs[33] = 1;
      buf[2].coeffs[32] = i;
      buf[2].coeffs[33] = 2;
      buf[3].coeffs[32] = i;
      buf[3].coeffs[33] = 3;
    }
    else {
      buf[0].coeffs[32] = 0;
      buf[0].coeffs[33] = i;
      buf[1].coeffs[32] = 1;
      buf[1].coeffs[33] = i;
      buf[2].coeffs[32] = 2;
      buf[2].coeffs[33] = i;
      buf[3].coeffs[32] = 3;
      buf[3].coeffs[33] = i;
    }

    shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

    ctr0 = rej_uniform_avx(a[i].vec[0].coeffs, buf[0].coeffs);
    ctr1 = rej_uniform_avx(a[i].vec[1].coeffs, buf[1].coeffs);
    ctr2 = rej_uniform_avx(a[i].vec[2].coeffs, buf[2].coeffs);
    ctr3 = rej_uniform_avx(a[i].vec[3].coeffs, buf[3].coeffs);

    while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
      shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

      ctr0 += rej_uniform(a[i].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
      ctr1 += rej_uniform(a[i].vec[1].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
      ctr2 += rej_uniform(a[i].vec[2].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
      ctr3 += rej_uniform(a[i].vec[3].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
    }

    poly_nttunpack(&a[i].vec[0]);
    poly_nttunpack(&a[i].vec[1]);
    poly_nttunpack(&a[i].vec[2]);
    poly_nttunpack(&a[i].vec[3]);
  }
}
#endif


#define gen_a(A,B)  gen_matrix(A,B,0)


//////////////////////////////////////////////













static void 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
uniform_vec_4x(int16_t *vec0, int16_t *vec1, int16_t *vec2, int16_t *vec3, size_t num, const uint8_t seed0[ZKPOP_SYMBYTES*2+5], const uint8_t seed1[ZKPOP_SYMBYTES*2+5], const uint8_t seed2[ZKPOP_SYMBYTES*2+5], const uint8_t seed3[ZKPOP_SYMBYTES*2+5])
{
  unsigned int ctr0, ctr1, ctr2, ctr3, buflen, off, k, j;
  uint8_t buf[4][SHAKE_RATE+2];
  keccakx4_state state;

  SHAKEx4_absorb_once(&state, seed0, seed1, seed2, seed3, ZKPOP_SYMBYTES*2+1);
  SHAKEx4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], 1, &state);

  buflen = SHAKE_RATE;
  ctr0 = rej_uniform_avx_len(vec0, num, buf[0], buflen);
  ctr1 = rej_uniform_avx_len(vec1, num, buf[1], buflen);
  ctr2 = rej_uniform_avx_len(vec2, num, buf[2], buflen);
  ctr3 = rej_uniform_avx_len(vec3, num, buf[3], buflen);

  while(ctr0 < num || ctr1 < num || ctr2 < num || ctr3 < num) 
  {
    off = buflen % 3;
    for (j = 0; j < 4; j++)
      for(k = 0; k < off; k++)
        buf[j][k] = buf[j][buflen - off + k];
    SHAKEx4_squeezeblocks(buf[0]+off, buf[1]+off, buf[2]+off, buf[3]+off, 1, &state);
    buflen = off + SHAKE_RATE;
    ctr0 += rej_uniform(vec0 + ctr0, num - ctr0, buf[0], buflen);
    ctr1 += rej_uniform(vec1 + ctr1, num - ctr1, buf[1], buflen);
    ctr2 += rej_uniform(vec2 + ctr2, num - ctr2, buf[2], buflen);
    ctr3 += rej_uniform(vec3 + ctr3, num - ctr3, buf[3], buflen);
  }
}

/*************************************************
* Name:        sample_audit_bundles
*
* Description: Sample ZKPOP_M-ZKPOP_SIGMA unique indices in [0, ZKPOP_M] that are
*              audited by the verifier. Emulated through ROM.
*
* Arguments:   - uint16_t *audit_bundles: output array (sorted indices)
*              - const uint8_t h[ZKPOP_SYMBYTES]: input value h (commitment to all values hij)
**************************************************/
#define ZKPOP_M_MASK (ZKPOP_M | (ZKPOP_M>>1) | (ZKPOP_M>>2) | (ZKPOP_M>>3) | (ZKPOP_M>>4) | (ZKPOP_M>>5) | (ZKPOP_M>>6) | (ZKPOP_M>>7) | (ZKPOP_M>>8) | (ZKPOP_M>>9) | (ZKPOP_M>>10) | (ZKPOP_M>>11) | (ZKPOP_M>>12) | (ZKPOP_M>>13) | (ZKPOP_M>>14) | (ZKPOP_M>>15))
void 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
sample_audit_bundles (uint8_t audit_bundles[(ZKPOP_M+7)/8], const uint8_t h[ZKPOP_SYMBYTES])
{
#if ZKPOP_M > 4096
#error "sample_audit_bundles is implemented for ZKPOP_M <= 2^12"
#endif
  size_t count = ZKPOP_SIGMA,i,offset = 0, len=0;
  keccak_state state;
  uint8_t buf[SHAKE_RATE+2];
  SHAKE_absorb_once(&state, h, ZKPOP_SYMBYTES);
  do {
    offset = len % 3;
    for (i = 0; i < offset; i++)
    {
      buf[i] = buf[len - offset + i];
    }
    SHAKE_squeezeblocks(buf + offset, 1, &state);
    len = offset + SHAKE_RATE;
    
    i = 0;
    // iterate over squeezed block (potentially drops some bits)
    while (len >= 3 && count < ZKPOP_M)
    {
      uint16_t tmp = (buf[i] | (buf[i+1] << 8)) & ZKPOP_M_MASK;
      if (tmp < count) // rejection sampling
      {
        audit_bundles[count/8] &= ~(1 << (count % 8));
        audit_bundles[count/8] |= ((audit_bundles[tmp / 8] >> (tmp % 8)) & 1) << (count % 8);
        audit_bundles[tmp / 8] |= 1 << (tmp % 8);
        
        count += 1;
      }
      if (count < ZKPOP_M)
      {
        tmp = ((buf[i+1] | (buf[i+2] << 8)) >> 4) & ZKPOP_M_MASK;
        if (tmp < count) // rejection sampling
        {
          audit_bundles[count/8] &= ~(1 << (count % 8));
          audit_bundles[count/8] |= (audit_bundles[tmp / 8] >> (tmp % 8)) << (count % 8);
          audit_bundles[tmp / 8] |= 1 << (tmp % 8);
          
          count += 1;
        }
      }
      len -= 3;
      i += 3;
    }
  } while(count < ZKPOP_M);
}

/*************************************************
* Name:        sample_audit_parties
*
* Description: Sample ZKPOP_TAU uniform random parties (range 0 to ZKPOP_N-1)
**************************************************/
void 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
sample_audit_parties (party_t audit_parties[ZKPOP_TAU], const uint8_t h2[ZKPOP_SYMBYTES])
{
  keccak_state state;
  SHAKE_absorb_once(&state, h2, ZKPOP_SYMBYTES);
/* ZKPOP_N might be a power of two (=no rejection sampling)
 * or not (=rejection sampling)
 */
#if ((ZKPOP_N-1)&ZKPOP_N) == 0 && ZKPOP_N <= 65536
  SHAKE_squeeze(audit_parties, ZKPOP_TAU*sizeof(party_t), &state);
  for (size_t j = 0; j < ZKPOP_TAU; j++) // ZKPOP_N is a power of two, so no rejection sampling is necessary
  {
    audit_parties[j] = audit_parties[j] % ZKPOP_N;
  }
#elif ((ZKPOP_N-1)&ZKPOP_N) != 0 && ZKPOP_N <= 256
#define ZKPOP_N_MASK (ZKPOP_N | (ZKPOP_N>>1) | (ZKPOP_N>>2) | (ZKPOP_N>>3) | (ZKPOP_N>>4) | (ZKPOP_N>>5) | (ZKPOP_N>>6) | (ZKPOP_N>>7))
  uint8_t shake_block[SHAKE_RATE];
  size_t num = 0,i;
  do {
    SHAKE_squeezeblocks(shake_block, 1, &state);
    for (i = 0; i < SHAKE_RATE; i++)
    {
      shake_block[i] &= ZKPOP_N_MASK;
      if (shake_block[i] < ZKPOP_N) // rejection sampling
      {
        audit_parties[num++] = shake_block[i];
        if (num >= ZKPOP_TAU)
        {
          break;
        }
      }
    }
  } while (num < ZKPOP_TAU);
#elif ((ZKPOP_N-1)&ZKPOP_N) != 0 && ZKPOP_N <= 65536
#define ZKPOP_N_MASK (ZKPOP_N | (ZKPOP_N>>1) | (ZKPOP_N>>2) | (ZKPOP_N>>3) | (ZKPOP_N>>4) | (ZKPOP_N>>5) | (ZKPOP_N>>6) | (ZKPOP_N>>7) | (ZKPOP_N>>8) | (ZKPOP_N>>9) | (ZKPOP_N>>10) | (ZKPOP_N>>11) | (ZKPOP_N>>12) | (ZKPOP_N>>13) | (ZKPOP_N>>14) | (ZKPOP_N>>15))
  uint8_t shake_block[SHAKE_RATE];
  size_t num = 0,i;
  do {
    SHAKE_squeezeblocks(shake_block, 1, &state);
    for (i = 0; i < SHAKE_RATE; i += 2)
    {
      audit_parties[num] = shake_block[i] | (((party_t)shake_block[i+1]) << 8);
      audit_parties[num] &= ZKPOP_N_MASK;
      if (audit_parties[num] < ZKPOP_N) // rejection sampling
      {
        num += 1;
        if (num >= ZKPOP_TAU)
        {
          break;
        }
      }
    }
  } while (num < ZKPOP_TAU);
#else
#error "Not yet implemented for the given choice of ZKPOP_N"
#endif
}

/*************************************************
* Name:        generate_seeds
*
* Description: Generate the seeds needed for generating the bundles vijk.
*
* Arguments:   - uint8_t *seedtree: pointer to the seed tree buffer
**************************************************/
void 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
generate_seeds(uint8_t *seedtree)
{
  size_t i, j, count;
  uint8_t *parent[4], *children[4];
  
  // expand initial seed to ZKPOP_TAU root seeds
  SHAKE(seedtree, ZKPOP_TAU*ZKPOP_SYMBYTES, seedtree, ZKPOP_SYMBYTES);
  
  // copy root seeds to correct position
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    memcpy(seedtree + (ZKPOP_TAU-j-1)*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES, seedtree + (ZKPOP_TAU-j-1)*ZKPOP_SYMBYTES, ZKPOP_SYMBYTES);
  }
  // build tree for each execution
  for (j = 0; j < ZKPOP_TAU; j += 4)
  {
    count = 1;
    for (i = 0; i < 4; i++)
    {
      parent[i] = &seedtree[((j+i) % ZKPOP_TAU)*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES];
      children[i] = &seedtree[((j+i) % ZKPOP_TAU)*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES + ZKPOP_SYMBYTES];
    }
    
    while (count < ZKPOP_N*2-1)
    {
      SHAKEx4(   children[0],
                 children[1],
                 children[2],
                 children[3],
                 ZKPOP_SYMBYTES * 2,
                 parent[0],
                 parent[1],
                 parent[2],
                 parent[3],
                 ZKPOP_SYMBYTES);
      for (i = 0; i < 4; i++)
      {
        children[i] += ZKPOP_SYMBYTES * 2;
        parent[i] += ZKPOP_SYMBYTES;
      }
      count += 2;
    }
  }
}

/*************************************************
* Name:        generate_S_E
*
* Description: Generate the polynomial vectors S and E from the small values vk that are not audited.
*
* Arguments:   - polyvec *S: output polynomial vector
*              - polyvec *E: output polynomial vector
*              - const bundle_t *vk: input bundles
*              - const uint8_t audit_bundles[(ZKPOP_M+7)/8]: indices of audited bundles (sorted)
**************************************************/
void 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
generate_S_E(polyvec *S, polyvec *E, const bundle_t *vk, const uint8_t audit_bundles[(ZKPOP_M+7)/8])
{
  size_t k, count = 0;
  uint8_t audit_flags = 0;
  for (k = 0; k < ZKPOP_M; k++)
  {
    if ((k%8) == 0)
    {
      audit_flags = audit_bundles[k / 8];
    }
    if (count < ZKPOP_M-ZKPOP_SIGMA && (audit_flags & 1))
    {
      count += 1;
    } else {
      if (k - count < KYBER_K*KYBER_N)
      {
        S->vec[(k-count)/KYBER_N].coeffs[(k-count)%KYBER_N] = vk->coeffs[k];
      } else {
        E->vec[(k-count-KYBER_K*KYBER_N)/KYBER_N].coeffs[(k-count-KYBER_K*KYBER_N)%KYBER_N] = vk->coeffs[k];
      }
    }
    audit_flags >>= 1;
  }
}

/*************************************************
* Name:        generate_S_E_openbundle
**************************************************/
void 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
generate_S_E_openbundle(polyvec *S, polyvec *E, open_bundle_t *openbundle, const bundle_t *vk, const uint8_t audit_bundles[(ZKPOP_M+7)/8])
{
  size_t k, count = 0;
  uint64_t audit_flags = 0;
  for (k = 0; k < ZKPOP_M; k++)
  {
    if ((k%64) == 0)
    {
      for (size_t i = 0; i < 8; i++)
      {
        audit_flags <<= 8;
        audit_flags |= audit_bundles[(k / 8) + (7-i)];
      }
    }
    if (count < ZKPOP_M-ZKPOP_SIGMA && (audit_flags & 1))
    {
      openbundle->coeffs[count] = vk->coeffs[k];
      count += 1;
    } else {
      if (k - count < KYBER_K*KYBER_N)
      {
        S->vec[(k-count)/KYBER_N].coeffs[(k-count)%KYBER_N] = vk->coeffs[k];
      } else {
        E->vec[(k-count-KYBER_K*KYBER_N)/KYBER_N].coeffs[(k-count-KYBER_K*KYBER_N)%KYBER_N] = vk->coeffs[k];
      }
    }
    audit_flags >>= 1;
  }
}

/*************************************************
* Name:        pack_nizkpop
*
* Description: Pack all values to the NIZKPoP buffer.
**************************************************/
static size_t 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
pack_nizkpop(uint8_t *zkpop, 
                           const uint8_t h1[ZKPOP_SYMBYTES], 
                           const uint8_t h2[ZKPOP_SYMBYTES], 
                           const uint8_t salt[ZKPOP_SYMBYTES], 
                           const uint8_t seedtree[(ZKPOP_N*2-1)*ZKPOP_TAU*ZKPOP_SYMBYTES], 
                           const uint8_t *comij,
                           const bundle_t deltabjk[ZKPOP_TAU],
                           const bundle_t *vk,
                           const party_t audit_parties[ZKPOP_TAU], 
                           const uint8_t audit_bundles[(ZKPOP_M+7)/8])
{
  size_t i, j, k, count=0;
  uint8_t audit_flags = 0;
  uint8_t *zkpop_cur = zkpop;
  size_t proof_size = KYBER_ZKPOP_MAXBYTES;
  
  STOREZKPOP(h1, ZKPOP_SYMBYTES);
  STOREZKPOP(h2, ZKPOP_SYMBYTES);
  STOREZKPOP(salt, ZKPOP_SYMBYTES);
  
  // store comij unaudited
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    STOREZKPOP(&comij[(audit_parties[j] * ZKPOP_TAU + j) * ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
  }
  
  // store all deltabjk
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (proof_size >= ZKPOP_M_PACKEDBYTES)
    {
      pack_bundle(zkpop_cur, &deltabjk[j]);
      zkpop_cur += ZKPOP_M_PACKEDBYTES;
      proof_size -= ZKPOP_M_PACKEDBYTES;
    } else {
      return __LINE__;
    }
  }
    
  // store all unaudited vk
  open_bundle_t openvk;
  for (k = 0; k < ZKPOP_M; k++)
  {
    if ((k%8) == 0)
    {
      audit_flags = audit_bundles[k/8];
    }
    if (audit_flags&1)
    {
      openvk.coeffs[count++] = vk->coeffs[k];
    } 
    audit_flags >>= 1;
  }
  for (k = 0; k < ZKPOP_M-ZKPOP_SIGMA; k += 16)
  {
    reduce16_avx(&openvk.vec[k/16], qdata.vec);
  }
  if (proof_size >= ZKPOP_M_SIGMA_PACKEDBYTES)
  {
    pack_open_bundle(zkpop_cur, &openvk);
    zkpop_cur += ZKPOP_M_SIGMA_PACKEDBYTES;
    proof_size -= ZKPOP_M_SIGMA_PACKEDBYTES;
  } else {
    return __LINE__;
  }
  
  // store seed tree nodes such that seedij[i][j] with i = audit_parties[j] can NOT be reconstructed
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    // send all seeds but one
    i = ZKPOP_N + audit_parties[j] - 1;
    while (i != 0)
    {
      // add sibling of current node
      if (i & 1)
      {
        STOREZKPOP(&seedtree[(j*(ZKPOP_N*2-1) + i+1)*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
      } else {
        STOREZKPOP(&seedtree[(j*(ZKPOP_N*2-1) + i-1)*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
      }
      // go to parent node
      i = (i-1) / 2;
    }
  }
  
  return KYBER_ZKPOP_MAXBYTES - proof_size;
}


/*************************************************
* Name:        build_subtree
*
* Description: Build subtree of a given node in the seed tree.
*
* Arguments:   - uint8_t *seedtree: pointer to the *BEGINNING* of the seed tree
*              - size_t start: offset of the parent node in the seed tree
**************************************************/
static void build_subtree(uint8_t *seedtree, size_t start) // TODO iterative and parallel
{
  if (2*start+1 < ZKPOP_N*2-1)
  {
    SHAKE(&seedtree[(2*start+1)*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES * 2, &seedtree[start*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
    build_subtree(seedtree, 2*start+1);
    build_subtree(seedtree, 2*start+2);
  }
}

/*************************************************
* Name:        unpack_nizkpop
*
* Description: Unpack given NIZKPoP buffer.
* 
* Returns 0 if ok.
**************************************************/
static int 
#ifdef PROFILE
__attribute__ ((noinline)) 
#endif
unpack_nizkpop(party_t audit_parties[ZKPOP_TAU],
                          uint8_t audit_bundles[(ZKPOP_M+7)/8],
                          open_bundle_t *vk,
                          uint8_t seedtree[(ZKPOP_N*2-1)*ZKPOP_TAU*ZKPOP_SYMBYTES],
                          bundle_t deltabjk[ZKPOP_TAU],
                          uint8_t salt[ZKPOP_SYMBYTES],
                          uint8_t h2[ZKPOP_SYMBYTES],
                          uint8_t h1[ZKPOP_SYMBYTES],
                          const uint8_t *zkpop,
                          size_t zkpop_size)
{
  size_t i,j;
  const uint8_t *zkpop_cur = zkpop;
  
  LOADZKPOP(h1, ZKPOP_SYMBYTES);
  LOADZKPOP(h2, ZKPOP_SYMBYTES);
  LOADZKPOP(salt, ZKPOP_SYMBYTES);
  
  // re-sample
  sample_audit_bundles(audit_bundles, h1);
  sample_audit_parties(audit_parties, h2);
  
  // unpack comij: we just skip this
  if (zkpop_size >= ZKPOP_TAU*ZKPOP_SYMBYTES)
  {
    zkpop_cur += ZKPOP_TAU*ZKPOP_SYMBYTES;
    zkpop_size -= ZKPOP_TAU*ZKPOP_SYMBYTES;
  } else {
    return -1;
  }
  
  // load all deltabjk
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (zkpop_size >= ZKPOP_M_PACKEDBYTES)
    {
      unpack_bundle(&deltabjk[j], zkpop_cur);
      zkpop_cur += ZKPOP_M_PACKEDBYTES;
      zkpop_size -= ZKPOP_M_PACKEDBYTES;
    } else {
      return __LINE__;
    }
  }
  
  // unpack all audited vk
  if (zkpop_size >= ZKPOP_M_SIGMA_PACKEDBYTES)
  {
    unpack_open_bundle(vk, zkpop_cur);
    zkpop_cur += ZKPOP_M_SIGMA_PACKEDBYTES;
    zkpop_size -= ZKPOP_M_SIGMA_PACKEDBYTES;
  } else {
    return __LINE__;
  }
  
  // unpack seed tree nodes such that seedij[i][j] where i = audit_parties[j] can NOT be reconstructed
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    i = ZKPOP_N + audit_parties[j] - 1;
    while (i != 0)
    {
      // load sibling of current node
      if (i & 1)
      {
        i += 1;
      } else {
        i -= 1;
      }
      LOADZKPOP(&seedtree[(j*(ZKPOP_N*2-1) + i)*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
      // re-build subtree
      build_subtree(&seedtree[(j*(ZKPOP_N*2-1) + 0)*ZKPOP_SYMBYTES], i);
      // go to parent node
      i = (i-1) / 2;
    }
  }
  
  return 0;
}



/*************************************************
* Name:        crypto_kem_keypair_nizkpop
*
* Description: Generate Kyber pk and sk as well a NIZKPoP of the sk.
*              Allocates memory for the NIZKPoP.
*
* Arguments:   - uint8_t pk[KYBER_PUBLICKEYBYTES]: output public key
*              - uint8_t sk[KYBER_SECRETKEYBYTES]: output secret key
*              - uint8_t **zkpop: output pointer to NIZKPoP buffer
*              - size_t *zkpop_size: size of NIZKPoP
* 
* Returns 0 if successful.
**************************************************/
int crypto_kem_keypair_nizkpop(uint8_t pk[KYBER_PUBLICKEYBYTES],
                       uint8_t sk[KYBER_SECRETKEYBYTES],
                       uint8_t **zkpop,
                       size_t *zkpop_size)
{
  // iterators
  size_t i, j, k, mi, ij;
  
  // ZKPOP values
  uint8_t salt[ZKPOP_SYMBYTES], h1[ZKPOP_SYMBYTES], h2[ZKPOP_SYMBYTES];
  uint8_t audit_bundles[(ZKPOP_M+7)/8] = {0};
  party_t audit_parties[ZKPOP_TAU];
#ifdef BIGPARAM
  uint8_t *seedtree = calloc(ZKPOP_TAU*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES, sizeof(uint8_t));
  if (seedtree == NULL)
  {
    return __LINE__;
  }
#else
  uint8_t seedtree[ZKPOP_TAU*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES];
#endif
#ifdef NO_RESAMPLING
  //bundle_t bijk[ZKPOP_N*ZKPOP_TAU];
  bundle_t *bijk = aligned_alloc(32, sizeof(bundle_t)*ZKPOP_N*ZKPOP_TAU);
  if (bijk == NULL)
  {
    return __LINE__;
  }
#else
  bundle_t bijk[4];
#endif
  bundle_t vk, deltabjk[ZKPOP_TAU] = {0};
#ifdef BIGPARAM
  uint8_t *comij = calloc(ZKPOP_N*ZKPOP_TAU*ZKPOP_SYMBYTES, sizeof(uint8_t));
  if (comij == NULL)
  {
    return __LINE__;
  }
#else
  uint8_t comij[ZKPOP_N*ZKPOP_TAU*ZKPOP_SYMBYTES];
#endif
  
  // ZKPoP buffers
  uint8_t comijbuf[4][ZKPOP_SYMBYTES*2 + 5];
  uint8_t packed_deltabjk[ZKPOP_M_PACKEDBYTES];
  uint8_t pvbuf[4][KYBER_POLYVECBYTES];
  uint8_t pvhash[4][ZKPOP_SYMBYTES];
  uint8_t openbundlebuf[4][ZKPOP_M_SIGMA_PACKEDBYTES];
  open_bundle_t openbundletmp;
  keccak_state state;
  
  // Kyber values
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  polyvec a[KYBER_K];
  polyvec B, S, E;
  
  
  // set proof size, allocate memory
  CALLOC8(*zkpop, KYBER_ZKPOP_MAXBYTES);
  
  // initial keygen randomness
  randombytes(buf, KYBER_SYMBYTES);
  hash_g(buf, buf, KYBER_SYMBYTES);
  randombytes(salt, ZKPOP_SYMBYTES);
  
  // sample vk from centered binomial distribution
  vec_cbd_eta1(&vk, noiseseed);
  
  // generate seed tree
  // root seeds for each tau generated from randomness
  randombytes(seedtree, ZKPOP_SYMBYTES);
  generate_seeds(seedtree);
  
  for (i = 0; i < 4; i++)
  {
    memcpy(comijbuf[i], salt, ZKPOP_SYMBYTES);
  }
  SHAKE_init(&state);
  SHAKE_absorb(&state, salt, ZKPOP_SYMBYTES);
  for (ij = 0; ij < ZKPOP_N*ZKPOP_TAU; ij += 4)
  {
    size_t nj = 4, I[4], J[4];
    if ((ij+4) >= ZKPOP_TAU*ZKPOP_N)
    {
      nj = ZKPOP_TAU*ZKPOP_N-ij;
    }
    // sample shares
    for (mi = 0; mi < 4; mi++)
    {
      I[mi] = (ij+mi)%ZKPOP_N;
      J[mi] = ((ij+mi)/ZKPOP_N)%ZKPOP_TAU;
      *((uint16_t*)(comijbuf[mi] + ZKPOP_SYMBYTES)) = I[mi];
      *((uint16_t*)(comijbuf[mi] + ZKPOP_SYMBYTES + 2)) = J[mi];
      comijbuf[mi][ZKPOP_SYMBYTES + 4] = 0;
      memcpy(comijbuf[mi] + ZKPOP_SYMBYTES + 5, &seedtree[(J[mi]*(ZKPOP_N*2-1) + ZKPOP_N + I[mi] - 1)*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
    }
#ifdef NO_RESAMPLING
    uniform_vec_4x( (int16_t*)bijk[J[0]*ZKPOP_N + I[0]].coeffs, 
                    (int16_t*)bijk[J[1]*ZKPOP_N + I[1]].coeffs, 
                    (int16_t*)bijk[J[2]*ZKPOP_N + I[2]].coeffs, 
                    (int16_t*)bijk[J[3]*ZKPOP_N + I[3]].coeffs, 
                    ZKPOP_M,
                    comijbuf[0],
                    comijbuf[1],
                    comijbuf[2],
                    comijbuf[3]);
#else
    uniform_vec_4x( (int16_t*)bijk[0].coeffs, 
                    (int16_t*)bijk[1].coeffs, 
                    (int16_t*)bijk[2].coeffs, 
                    (int16_t*)bijk[3].coeffs, 
                    ZKPOP_M,
                    comijbuf[0],
                    comijbuf[1],
                    comijbuf[2],
                    comijbuf[3]);
#endif
    
    // commit to seeds
    for (mi = 0; mi < 4; mi++)
    {
      comijbuf[mi][ZKPOP_SYMBYTES + 4] = 1;
    }
    SHAKEx4( 
      &comij[(I[0] * ZKPOP_TAU + J[0]) * ZKPOP_SYMBYTES],
      &comij[(I[1] * ZKPOP_TAU + J[1]) * ZKPOP_SYMBYTES],
      &comij[(I[2] * ZKPOP_TAU + J[2]) * ZKPOP_SYMBYTES],
      &comij[(I[3] * ZKPOP_TAU + J[3]) * ZKPOP_SYMBYTES],
      ZKPOP_SYMBYTES,
      comijbuf[0],
      comijbuf[1],
      comijbuf[2],
      comijbuf[3],
      ZKPOP_SYMBYTES*2+5);
    
    // update deltabjk
    for (mi = 0; mi < nj; mi++)
    {
      for (k = 0; k < ZKPOP_M; k += 16)
      {
#ifdef NO_RESAMPLING
        deltabjk[J[mi]].vec[k/16] = _mm256_add_epi16(deltabjk[J[mi]].vec[k/16], bijk[J[mi]*ZKPOP_N + I[mi]].vec[k/16]);
#else
        deltabjk[J[mi]].vec[k/16] = _mm256_add_epi16(deltabjk[J[mi]].vec[k/16], bijk[mi].vec[k/16]);
#endif
        reduce16_avx(&deltabjk[J[mi]].vec[k/16], qdata.vec);
      }
    }
    
    for (mi = 0; mi < nj; mi++)
    {
      SHAKE_absorb(&state, &comij[(I[mi] * ZKPOP_TAU + J[mi]) * ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
    }
  }
  
  // compute final deltabjk and finalize h1
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    for (k = 0; k < ZKPOP_M; k += 16)
    {
      deltabjk[j].vec[k/16] = _mm256_sub_epi16(vk.vec[k/16], deltabjk[j].vec[k/16]);
      reduce16_avx(&deltabjk[j].vec[k/16], qdata.vec);
    }
    pack_bundle(packed_deltabjk, &deltabjk[j]);
    SHAKE_absorb(&state, packed_deltabjk, ZKPOP_M_PACKEDBYTES);
  }
  
  // here, we could additionally absorb attributes
  SHAKE_finalize(&state);
  SHAKE_squeeze(h1, ZKPOP_SYMBYTES, &state);
  
  // sample audit_bundles
  sample_audit_bundles(audit_bundles, h1);
  
  // generate matrix A as in original key gen
  gen_a(a, publicseed);
  
  // generate Sij, Eij from those vijk values that are selected for auditing
  SHAKE_init(&state);
  SHAKE_absorb(&state, salt, ZKPOP_SYMBYTES);
  SHAKE_absorb(&state, h1, ZKPOP_SYMBYTES);
  keccak_state bundlestate;
  SHAKE_init(&bundlestate);
  
  // compute Bij = A*Sij + Eij
  for (ij = 0; ij < ZKPOP_N*ZKPOP_TAU; ij += 4)
  {
    // re-sample bijk
    size_t I[4], J[4], ijcnt = 4;
    
    if ((ij+4) >= ZKPOP_N*ZKPOP_TAU)
    {
      ijcnt = ZKPOP_N*ZKPOP_TAU - ij;
    }
    
#ifdef NO_RESAMPLING
    for (mi = 0; mi < ijcnt; mi++)
    {
      I[mi] = (ij+mi)%ZKPOP_N;
      J[mi] = (ij+mi)/ZKPOP_N;
    }
#else
    // sample shares
    for (mi = 0; mi < ijcnt; mi++)
    {
      I[mi] = (ij+mi)%ZKPOP_N;
      J[mi] = (ij+mi)/ZKPOP_N;
      *((uint16_t*)(comijbuf[mi] + ZKPOP_SYMBYTES)) = I[mi];
      *((uint16_t*)(comijbuf[mi] + ZKPOP_SYMBYTES + 2)) = J[mi];
      comijbuf[mi][ZKPOP_SYMBYTES + 4] = 0;
      memcpy(comijbuf[mi] + ZKPOP_SYMBYTES + 5, &seedtree[(J[mi]*(ZKPOP_N*2-1) + ZKPOP_N + I[mi] - 1)*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
    }
    uniform_vec_4x( bijk[0].coeffs, 
                    bijk[1].coeffs, 
                    bijk[2].coeffs, 
                    bijk[3].coeffs, 
                    ZKPOP_M,
                    comijbuf[0],
                    comijbuf[1],
                    comijbuf[2],
                    comijbuf[3]);
#endif
    
    for (mi = 0; mi < ijcnt; mi++)
    {
      // update first party
      if (I[mi] == 0)
      {
        for (k = 0; k < ZKPOP_M; k += 16)
        {
#ifdef NO_RESAMPLING
          bijk[J[mi]*ZKPOP_N + I[mi]].vec[k/16] = _mm256_add_epi16(bijk[J[mi]*ZKPOP_N + I[mi]].vec[k/16], deltabjk[J[mi]].vec[k/16]);
          reduce16_avx(&bijk[J[mi]*ZKPOP_N + I[mi]].vec[k/16], qdata.vec);
#else
          bijk[mi].vec[k/16] = _mm256_add_epi16(bijk[mi].vec[k/16], deltabjk[J[mi]].vec[k/16]);
          reduce16_avx(&bijk[mi].vec[k/16], qdata.vec);
#endif
        }
      }
      
      // generate Sij and Eij from bundle
#ifdef NO_RESAMPLING
      generate_S_E_openbundle(&S, &E, &openbundletmp, &bijk[J[mi]*ZKPOP_N + I[mi]], audit_bundles);
#else
      generate_S_E_openbundle(&S, &E, &openbundletmp, &bijk[mi], audit_bundles);
#endif
      pack_open_bundle(openbundlebuf[mi], &openbundletmp);
      
      // compute Bij = A*Sij + Eij
      polyvec_ntt(&S);
      polyvec_reduce(&S);
      polyvec_ntt(&E);
      for(k = 0; k < KYBER_K; k++) 
      {
        polyvec_basemul_acc_montgomery(&B.vec[k], &a[k], &S);
        poly_tomont(&B.vec[k]);
      }

      polyvec_add(&B, &B, &E);
      polyvec_reduce(&B);
    
      // hash result
      polyvec_tobytes(pvbuf[mi], &B);
    }
    SHAKEx4(
      pvhash[0], 
      pvhash[1],
      pvhash[2],
      pvhash[3],
      ZKPOP_SYMBYTES,
      pvbuf[0],
      pvbuf[1],
      pvbuf[2],
      pvbuf[3],
      KYBER_POLYVECBYTES);
    SHAKE_absorb(&state, pvhash[0], ZKPOP_SYMBYTES*ijcnt);
    
    SHAKEx4(
      pvhash[0], 
      pvhash[1],
      pvhash[2],
      pvhash[3],
      ZKPOP_SYMBYTES,
      openbundlebuf[0],
      openbundlebuf[1],
      openbundlebuf[2],
      openbundlebuf[3],
      ZKPOP_M_SIGMA_PACKEDBYTES);
    SHAKE_absorb(&bundlestate, pvhash[0], ZKPOP_SYMBYTES*ijcnt);
  }
  
  // generate S and E from those vk values that are NOT selected for auditing
  generate_S_E(&S, &E, &vk, audit_bundles);
  
  // compute B as in original key gen
  polyvec_ntt(&S);
  polyvec_reduce(&S);
  polyvec_ntt(&E);
  for(i = 0; i < KYBER_K; i++) 
  {
    polyvec_basemul_acc_montgomery(&B.vec[i], &a[i], &S);
    poly_tomont(&B.vec[i]);
  }

  polyvec_add(&B, &B, &E);
  polyvec_reduce(&B);
  
  // hash result
  polyvec_tobytes(pvbuf[0], &B);
  SHAKE_absorb(&state, pvbuf[0], KYBER_POLYVECBYTES);
  SHAKE_finalize(&bundlestate);
  SHAKE_squeeze(h2, ZKPOP_SYMBYTES, &bundlestate); // use h2 to finalize the audited bundles hash
  SHAKE_absorb(&state, h2, ZKPOP_SYMBYTES);
  SHAKE_finalize(&state);
  SHAKE_squeeze(h2, ZKPOP_SYMBYTES, &state);
  
  // sample audit_parties
  sample_audit_parties(audit_parties, h2);
  
  // pack the NIZKPoP
  *zkpop_size = pack_nizkpop(*zkpop, h1, h2, salt, seedtree, comij, deltabjk, &vk, audit_parties, audit_bundles);
  
#ifdef BIGPARAM
  free(comij);
  free(seedtree);
#endif
#ifdef NO_RESAMPLING
  free(bijk);
#endif
  
  pack_sk(sk, &S); // S
  pack_pk(pk, &B, publicseed); // B
  
  memcpy(sk+KYBER_INDCPA_SECRETKEYBYTES, pk, KYBER_INDCPA_PUBLICKEYBYTES);
  hash_h(sk+KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES);
  /* Value z for pseudo-random output on reject */
  randombytes(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES);
  return 0;
}

/*************************************************
* Name:        crypto_nizkpop_verify
*
* Description: Verify a NIZKPoP for a given pk.
*
* Arguments:   - const unsigned char pk[KYBER_PUBLICKEYBYTES]: public key
*              - const unsigned char *zkpop: NIZKPoP buffer
*              - unsigned long zkpop_size: NIZKPoP size
**************************************************/
int crypto_nizkpop_verify(const unsigned char pk[KYBER_PUBLICKEYBYTES], const unsigned char *zkpop, unsigned long zkpop_size)
{
  // iterators
  size_t i,j,k,mi,i2;
  size_t count = 0, I[4] = {0};
  
  // ZKPoP values and buffers
  party_t audit_parties[ZKPOP_TAU] = {0};
  uint8_t audit_bundles[(ZKPOP_M+7)/8] = {0};
  uint8_t h1[ZKPOP_SYMBYTES], h2[ZKPOP_SYMBYTES], salt[ZKPOP_SYMBYTES], h1_check[ZKPOP_SYMBYTES], h2_check[ZKPOP_SYMBYTES];
#ifndef BIGPARAM
  uint8_t seedtree[ZKPOP_TAU*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES] = {0};
#else
  uint8_t *seedtree = calloc(ZKPOP_TAU*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES, sizeof(uint8_t));
#endif
  open_bundle_t vk, missingshare, tmpshare;
  bundle_t bijk[4] = {0}, deltabjk[ZKPOP_TAU];
  polyvec Eij;
  polyvec Sij;
  polyvec Bij, missingBij;
  uint8_t Bijhash[ZKPOP_N][ZKPOP_SYMBYTES], openbundlehash[ZKPOP_N][ZKPOP_SYMBYTES], tmphash[4][ZKPOP_SYMBYTES];
  uint8_t pvbuf[4][KYBER_POLYVECBYTES];
  uint8_t comijbuf[4][ZKPOP_SYMBYTES*2 + 5] = {0};
  uint8_t openbundlebuf[4][ZKPOP_M_SIGMA_PACKEDBYTES];
  keccak_state h1state, h2state, bundlestate;
  const uint8_t *zkpptr = zkpop + 3*ZKPOP_SYMBYTES; // used for absorbing the transmitted comij
  uint8_t packed_deltabjk[ZKPOP_M_PACKEDBYTES];
  
  // Kyber values
  uint8_t publicseed[KYBER_SYMBYTES];
  polyvec B;
  polyvec a[KYBER_K];
  
  unpack_pk(&B, publicseed, pk);
  
  
  if (unpack_nizkpop( audit_parties,
                      audit_bundles,
                      &vk,
                      seedtree,
                      deltabjk,
                      salt,
                      h2, 
                      h1, 
                      zkpop, 
                      zkpop_size
                     ) != 0)
  {
    return __LINE__;
  }
  
  // generate matrix A
  gen_a(a, publicseed);
  
  // re-compute h1 and h2
  SHAKE_init(&bundlestate);
  SHAKE_init(&h1state);
  SHAKE_init(&h2state);
  SHAKE_absorb(&h2state, salt, ZKPOP_SYMBYTES);
  SHAKE_absorb(&h1state, salt, ZKPOP_SYMBYTES);
  SHAKE_absorb(&h2state, h1, ZKPOP_SYMBYTES);
  
  for (i = 0; i < 4; i++)
  {
    memcpy(comijbuf[i], salt, ZKPOP_SYMBYTES);
  }
  
  // compute Bij = A*Sij + Eij
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    memcpy(&missingBij, &B, sizeof(polyvec));
  
    memcpy(&missingshare, &vk, sizeof(open_bundle_t));
    
    for (i = 0; i < ZKPOP_N; i++)
    {
      if (i == audit_parties[j])
      {
        if (count == 0)
        {
          SHAKE_absorb(&h1state, zkpptr, ZKPOP_SYMBYTES);
          zkpptr += ZKPOP_SYMBYTES;
        }
        continue;
      }
      I[count] = i;
      *((uint16_t*)(comijbuf[count] + ZKPOP_SYMBYTES)) = I[count];
      *((uint16_t*)(comijbuf[count] + ZKPOP_SYMBYTES + 2)) = j;
      comijbuf[count][ZKPOP_SYMBYTES + 4] = 0;
      memcpy(comijbuf[count] + ZKPOP_SYMBYTES + 5, &seedtree[(j*(ZKPOP_N*2-1) + ZKPOP_N + I[count] - 1)*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
      count += 1;
      
      if (count == 4 || i == ZKPOP_N-1 || (i == ZKPOP_N-2 && ZKPOP_N-1 == audit_parties[j]))
      {
        uniform_vec_4x( bijk[0].coeffs, 
                        bijk[1].coeffs, 
                        bijk[2].coeffs, 
                        bijk[3].coeffs, 
                        ZKPOP_M,
                        comijbuf[0],
                        comijbuf[1],
                        comijbuf[2],
                        comijbuf[3]);
        
        for (mi = 0; mi < count; mi++)
        {
          // update first party
          if (I[mi] == 0)
          {
            for (k = 0; k < ZKPOP_M; k += 16)
            {
              bijk[mi].vec[k/16] = _mm256_add_epi16(bijk[mi].vec[k/16], deltabjk[j].vec[k/16]);
              reduce16_avx(&bijk[mi].vec[k/16], qdata.vec);
            }
          }
          
          // generate Sij and Eij from bundle
          generate_S_E_openbundle(&Sij, &Eij, &tmpshare, &bijk[mi], audit_bundles);
          
          // compute Bij = A*Sij + Eij
          polyvec_ntt(&Sij);
          polyvec_reduce(&Sij);
          polyvec_ntt(&Eij);
          for(k = 0; k < KYBER_K; k++) 
          {
            polyvec_basemul_acc_montgomery(&Bij.vec[k], &a[k], &Sij);
            poly_tomont(&Bij.vec[k]);
          }

          polyvec_add(&Bij, &Bij, &Eij);
          polyvec_reduce(&Bij);
          polyvec_tobytes(pvbuf[mi], &Bij);
          
          for (k = 0; k < KYBER_K; k++)
          {
            poly_sub(&missingBij.vec[k], &missingBij.vec[k], &Bij.vec[k]);
          }
          polyvec_reduce(&missingBij);
          
          // update missingshare
          for (k = 0; k < ZKPOP_M-ZKPOP_SIGMA; k += 16)
          {
            missingshare.vec[k/16] = _mm256_sub_epi16(missingshare.vec[k/16], tmpshare.vec[k/16]);
            reduce16_avx(&missingshare.vec[k/16], qdata.vec);
          }
          pack_open_bundle(openbundlebuf[mi], &tmpshare);
        }
        
        SHAKEx4(
          tmphash[0],
          tmphash[1],
          tmphash[2],
          tmphash[3],
          ZKPOP_SYMBYTES,
          pvbuf[0],
          pvbuf[1],
          pvbuf[2],
          pvbuf[3],
          KYBER_POLYVECBYTES
        );
        
        for (mi = 0; mi < count; mi++)
        {
          memcpy(Bijhash[I[mi]], tmphash[mi], ZKPOP_SYMBYTES);
        }
        
        SHAKEx4(
          tmphash[0],
          tmphash[1],
          tmphash[2],
          tmphash[3],
          ZKPOP_SYMBYTES,
          openbundlebuf[0],
          openbundlebuf[1],
          openbundlebuf[2],
          openbundlebuf[3],
          ZKPOP_M_SIGMA_PACKEDBYTES
        );
        
        for (mi = 0; mi < count; mi++)
        {
          memcpy(openbundlehash[I[mi]], tmphash[mi], ZKPOP_SYMBYTES);
        }
        
        // update h1
        for (mi = 0; mi < count; mi++)
        {
          comijbuf[mi][ZKPOP_SYMBYTES + 4] = 1;
        }
        SHAKEx4(
          tmphash[0],
          tmphash[1],
          tmphash[2],
          tmphash[3],
          ZKPOP_SYMBYTES,
          comijbuf[0],
          comijbuf[1],
          comijbuf[2],
          comijbuf[3],
          ZKPOP_SYMBYTES*2+5
        );
        for (i2 = I[0],mi = 0; i2 <= I[count-1]; i2++)
        {
          if (i2 == audit_parties[j])
          {
            SHAKE_absorb(&h1state, zkpptr, ZKPOP_SYMBYTES);
            zkpptr += ZKPOP_SYMBYTES;
          } else {
            SHAKE_absorb(&h1state, tmphash[mi++], ZKPOP_SYMBYTES);
          }
        }
        
        count = 0;
      }
    }
    
    polyvec_tobytes(pvbuf[0], &missingBij);
    SHAKE(Bijhash[audit_parties[j]], ZKPOP_SYMBYTES, pvbuf[0], KYBER_POLYVECBYTES);
    SHAKE_absorb(&h2state, Bijhash[0], ZKPOP_SYMBYTES*ZKPOP_N);
    
    pack_open_bundle(openbundlebuf[0], &missingshare);
    SHAKE(openbundlehash[audit_parties[j]], ZKPOP_SYMBYTES, openbundlebuf[0], ZKPOP_M_SIGMA_PACKEDBYTES);
    SHAKE_absorb(&bundlestate, openbundlehash[0], ZKPOP_SYMBYTES*ZKPOP_N);
  }
  
  polyvec_tobytes(pvbuf[0], &B);
  SHAKE_absorb(&h2state, pvbuf[0], KYBER_POLYVECBYTES);
  
  SHAKE_finalize(&bundlestate);
  SHAKE_squeeze(tmphash[0], ZKPOP_SYMBYTES, &bundlestate);
  SHAKE_absorb(&h2state, tmphash[0], ZKPOP_SYMBYTES);
  SHAKE_finalize(&h2state);
  SHAKE_squeeze(h2_check, ZKPOP_SYMBYTES, &h2state);
  
  if (memcmp(h2, h2_check, ZKPOP_SYMBYTES) != 0)
  {
    return __LINE__;
  }
  
  // finalize h1
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    pack_bundle(packed_deltabjk, &deltabjk[j]);
    SHAKE_absorb(&h1state, packed_deltabjk, ZKPOP_M_PACKEDBYTES);
  }
  SHAKE_finalize(&h1state);
  SHAKE_squeeze(h1_check, ZKPOP_SYMBYTES, &h1state);
  
  if (memcmp(h1, h1_check, ZKPOP_SYMBYTES) != 0)
  {
    return __LINE__;
  }
  
  
#ifdef BIGPARAM
  free(seedtree);
#endif
  
  return 0;
}
