#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "symmetric.h"
#include "randombytes.h"
#include "cbd.h"
#include "reduce.h"
#include <stdio.h>
#include "fips202.h"

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *r: pointer to the output serialized public key
*              polyvec *pk: pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  size_t i;
  polyvec_tobytes(r, pk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    r[i+KYBER_POLYVECBYTES] = seed[i];
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
  size_t i;
  polyvec_frombytes(pk, packedpk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    seed[i] = packedpk[i+KYBER_POLYVECBYTES];
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *r: pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key; inverse of pack_sk
*
* Arguments:   - polyvec *sk: pointer to output vector of polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk(polyvec *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        pack_ciphertext
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk: pointer to the input vector of polynomials b
*              poly *v: pointer to the input polynomial v
**************************************************/
static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES], polyvec *b, poly *v)
{
  polyvec_compress(r, b);
  poly_compress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b: pointer to the output vector of polynomials b
*              - poly *v: pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext(polyvec *b, poly *v, const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  poly_decompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output buffer
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
  while(ctr < len && pos + 3 <= buflen) {
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

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a: pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed: boolean deciding whether A or A^T is generated
**************************************************/
#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
// Not static for benchmarking
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j, k;
  unsigned int buflen, off;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+2];
  xof_state state;

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;
      ctr = rej_uniform(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

      while(ctr < KYBER_N) {
        off = buflen % 3;
        for(k = 0; k < off; k++)
          buf[k] = buf[buflen - off + k];
        xof_squeezeblocks(buf + off, 1, &state);
        buflen = off + XOF_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}

/*************************************************
* Name:        uniform_vec
*
* Description: Sample num uniform coefficients mod q
*
* Arguments:   - int16_t *vec: pointer to ouptput array
*              - size_t num: number of coefficients to sample
*              - const uint8_t *seed: pointer to input seed
**************************************************/
static void uniform_vec(int16_t *vec, size_t num, const uint8_t seed[32])
{
  unsigned int ctr, k;
  unsigned int buflen, off;
  uint8_t buf[XOF_BLOCKBYTES+2];
  xof_state state;

  xof_absorb(&state, seed, 0, 0);
  
  xof_squeezeblocks(buf, 1, &state);
  buflen = XOF_BLOCKBYTES;
  ctr = rej_uniform(vec, num, buf, buflen);

  while(num > ctr) 
  {
    off = buflen % 3;
    for(k = 0; k < off; k++)
      buf[k] = buf[buflen - off + k];
    xof_squeezeblocks(buf + off, 1, &state);
    buflen = off + XOF_BLOCKBYTES;
    ctr += rej_uniform(vec + ctr, num - ctr, buf, buflen);
  }
}

// TODO are there any more efficient options?
#define ZKPOP_M_MASK (ZKPOP_M | (ZKPOP_M>>1) | (ZKPOP_M>>2) | (ZKPOP_M>>3) | (ZKPOP_M>>4) | (ZKPOP_M>>5) | (ZKPOP_M>>6) | (ZKPOP_M>>7) | (ZKPOP_M>>8) | (ZKPOP_M>>9) | (ZKPOP_M>>10) | (ZKPOP_M>>11) | (ZKPOP_M>>12) | (ZKPOP_M>>13) | (ZKPOP_M>>14) | (ZKPOP_M>>15))
static void sample_audit_bundles (uint16_t *audit_bundles, const uint8_t h[32])
{
  size_t count = 0,i,k;
  keccak_state state;
  uint8_t buf[SHAKE256_RATE];
  shake256_absorb_once(&state, h, 32);
  do {
    shake256_squeezeblocks(buf, 1, &state);
    
    // iterate over squeezed block (this drops some bits)
    for (i = 0; i < SHAKE256_RATE/2 && count < ZKPOP_M-ZKPOP_SIGMA; i++)
    {
      uint16_t tmp = *((uint16_t*)&buf[2*i]) & ZKPOP_M_MASK;
      if (tmp < ZKPOP_M)
      {
        if (count > 0)
        {
          long long int index_l = 0;
          long long int index_h = (long long int)count-1;
          while (index_l < index_h) // binary search
          {
            size_t mid = (index_l + index_h) / 2;
            if (audit_bundles[mid] < tmp)
            {
              index_l = (long long int)mid + 1;
            } else if (audit_bundles[mid] == tmp) {
              index_l = mid;
              break;
            } else {
              index_h = (long long int)mid - 1;
            }
          }
          if (audit_bundles[index_l] != tmp) // insert if not found
          {
            for (k = count - index_l; k > 0; k--)
            {
              audit_bundles[index_l + k] = audit_bundles[index_l + k - 1];
            }
            if (audit_bundles[index_l] < tmp)
            {
              audit_bundles[index_l+1] = tmp;
            } else {
              audit_bundles[index_l] = tmp;
            }
            count += 1;
          }
        } else {
          audit_bundles[0] = tmp;
          count = 1;
        }
      }
    }
  } while(count < ZKPOP_M-ZKPOP_SIGMA);
}

static void sample_audit_parties (uint16_t audit_parties[ZKPOP_TAU], const uint8_t *h, const uint8_t *hB, const uint8_t publicseed[KYBER_SYMBYTES], const polyvec *B_true)
{
  uint8_t buf[32 + 32 + KYBER_SYMBYTES + sizeof(polyvec)];
  memcpy(buf, h, 32);
  memcpy(buf + 32, hB, 32);
  memcpy(buf + 64, publicseed, KYBER_SYMBYTES);
  memcpy(buf + 64 + KYBER_SYMBYTES, B_true, sizeof(polyvec));
/* ZKPOP_N might be a power of two (=no rejection sampling)
 * or not (=rejection sampling)
 * Also, ZKPOP_N might be <= 8 bit (working on bytes squeezed out of shake)
 * or >8 bit (working on 16-bit words squeezed out of shake)
 */
#if ((ZKPOP_N-1)&ZKPOP_N) == 0 && ZKPOP_N <= 256
  uint8_t tmp[ZKPOP_TAU];
  shake256(tmp, ZKPOP_TAU, buf, 32 + 32 + KYBER_SYMBYTES + sizeof(polyvec));
  for (size_t j = 0; j < ZKPOP_TAU; j++) // ZKPOP_N is a power of two, so no rejection sampling is necessary
  {
    audit_parties[j] = tmp[j] % ZKPOP_N;
  }
#elif ((ZKPOP_N-1)&ZKPOP_N) == 0 && ZKPOP_N <= 65536
  shake256((uint8_t*)audit_parties, ZKPOP_TAU * sizeof(uint16_t), buf, 32 + 32 + KYBER_SYMBYTES + sizeof(polyvec));
  for (size_t j = 0; j < ZKPOP_TAU; j++) // ZKPOP_N is a power of two, so no rejection sampling is necessary
  {
    audit_parties[j] %= ZKPOP_N;
  }
#elif ((ZKPOP_N-1)&ZKPOP_N) != 0 && ZKPOP_N <= 256
#define ZKPOP_N_MASK (ZKPOP_N | (ZKPOP_N>>1) | (ZKPOP_N>>2) | (ZKPOP_N>>3) | (ZKPOP_N>>4) | (ZKPOP_N>>5) | (ZKPOP_N>>6) | (ZKPOP_N>>7))
  uint8_t shake_block[SHAKE256_RATE];
  size_t num = 0,i;
  keccak_state state;
  shake256_init(&state);
  shake256_absorb(&state, buf, 32 + 32 + KYBER_SYMBYTES + sizeof(polyvec));
  shake256_finalize(&state);
  do {
    shake256_squeezeblocks(shake_block, 1, &state);
    for (i = 0; i < SHAKE256_RATE/8; i++)
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
  uint8_t shake_block[SHAKE256_RATE];
  size_t num = 0,i;
  keccak_state state;
  shake256_init(&state);
  shake256_absorb(&state, buf, 32 + 32 + KYBER_SYMBYTES + sizeof(polyvec));
  shake256_finalize(&state);
  do {
    shake256_squeezeblocks(shake_block, 1, &state);
    for (i = 0; i < SHAKE256_RATE/16; i++)
    {
      uint16_t *tmp = (uint16_t*)(&shake_block[2*i]);
      *tmp &= ZKPOP_N_MASK;
      if (*tmp < ZKPOP_N) // rejection sampling
      {
        audit_parties[num++] = *tmp;
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

#define LOADZKPOP(dest, size) do {if (zkpop_size < size) { printf("bad zkpop: size\n"); return __LINE__;} memcpy(dest, zkpop_cur, size); zkpop_cur += size; zkpop_size -= size;} while(0)
#define LOADZKPOP_POLYVEC(dest) do {if (zkpop_size < KYBER_POLYVECBYTES) { printf("bad zkpop: size\n"); return __LINE__;} polyvec_frombytes(dest, zkpop_cur); zkpop_cur += KYBER_POLYVECBYTES; zkpop_size -= KYBER_POLYVECBYTES;} while(0)
#define STOREZKPOP(src, size) do {if (proof_size < size) { printf("bad zkpop: size\n"); return __LINE__;} memcpy(zkpop_cur, src, size); zkpop_cur += size; proof_size -= size;} while(0)
#define STOREZKPOP_POLYVEC(src) do {if (proof_size < KYBER_POLYVECBYTES) { printf("bad zkpop: size\n"); return __LINE__;} polyvec_tobytes(zkpop_cur, src); zkpop_cur += KYBER_POLYVECBYTES; proof_size -= KYBER_POLYVECBYTES;} while(0)


#define CALLOCPOLYVEC(x, size) do { x = calloc((size_t)size, sizeof(polyvec)); if (x == NULL) {printf("bad zkpop: bad alloc\n");return __LINE__;} } while(0)
#define CALLOC16(x, size) do { x = calloc((size_t)size, sizeof(uint16_t)); if (x == NULL) {printf("bad zkpop: bad alloc\n");return __LINE__;} } while(0)
#define CALLOC8(x, size) do { x = calloc((size_t)size, sizeof(uint8_t)); if (x == NULL) {printf("bad zkpop: bad alloc\n");return __LINE__;} } while(0)
/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
int indcpa_keypair_zkpop(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES],
                    uint8_t **zkpop,
                    size_t *zkpop_size)
{
  size_t i, j, k, count, mi, proof_size;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  
  uint16_t audit_bundles[ZKPOP_M-ZKPOP_SIGMA];
  uint8_t rootseeds[ZKPOP_TAU][32];
  uint8_t *seedtree;
  int16_t *vijk;
  uint8_t *hij, *hBij;
  uint8_t h[32], hB[32];
  polyvec *B, *S, *E;
  polyvec B_true, S_true, E_true;
  uint16_t audit_parties[ZKPOP_TAU];
  
/**
 * For KYBER_ETA1 == 2, min(max_len,8*(len/4)) samples are written to vec.
 * For KYBER_ETA1 == 3, min(max_len,4*(len/3)) elements are written to vec.
 */
#if KYBER_ETA1 == 2
  uint8_t samplebuf[(ZKPOP_M*4+7)/8];
#elif KYBER_ETA1 == 3
  uint8_t samplebuf[(ZKPOP_M*3+3)/4];
#endif
  int16_t vk[ZKPOP_M];
  
  polyvec a[KYBER_K];

  randombytes(buf, KYBER_SYMBYTES);
  hash_g(buf, buf, KYBER_SYMBYTES);
  
  // sample vk
#if KYBER_ETA1 == 2
  prf(samplebuf, (ZKPOP_M*4+7)/8, noiseseed, KYBER_SYMBYTES);
  vec_cbd_eta1(vk, ZKPOP_M, samplebuf, (ZKPOP_M*4+7)/8);
#elif KYBER_ETA1 == 3
  prf(samplebuf, (ZKPOP_M*3+3)/4, noiseseed, KYBER_SYMBYTES);
  vec_cbd_eta1(vk, ZKPOP_M, samplebuf, (ZKPOP_M*3+3)/4);
#endif
  
  
  // generate seed tree
  CALLOC8(seedtree, ZKPOP_TAU*((ZKPOP_N-1)*2-1)*32);
  // root seeds for each tau generated from randomness
  randombytes(rootseeds[0], 32);
  shake256(rootseeds[0], ZKPOP_TAU*32, rootseeds[0], 32);
  // build tree for each execution
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    count = 1;
    uint8_t *parent, *children;
    memcpy(&seedtree[j*((ZKPOP_N-1)*2-1)*32], rootseeds[j], 32);
    parent = &seedtree[j*((ZKPOP_N-1)*2-1)*32];
    children = &seedtree[j*((ZKPOP_N-1)*2-1)*32 + 32];
    
    while (count < (ZKPOP_N-1)*2-1)
    {
      prf(children, 64, parent, 0); // either SHAKE256 or AES256
      children += 64;
      parent += 32;
      count += 2;
    }
  }
  
  // sample vijk
  CALLOC16(vijk, ZKPOP_N*ZKPOP_TAU*ZKPOP_M);
  for (i = 0; i < ZKPOP_N-1; i++)
  {
    for (j = 0; j < ZKPOP_TAU; j++)
    {
      uniform_vec(&vijk[(i*ZKPOP_TAU + j) * ZKPOP_M], ZKPOP_M, &seedtree[(j*((ZKPOP_N-1)*2-1) + ZKPOP_N + i - 2)*32]);
    }
  }
  
  // compute vNjk
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    for (k = 0; k < ZKPOP_M; k++)
    {
      int16_t tmp = 0;
      for (i = 0; i < ZKPOP_N-1; i++)
      {
        tmp += vijk[(i*ZKPOP_TAU + j) * ZKPOP_M + k];
        if (((i+1)%9) == 0)
        {
          tmp = barrett_reduce(tmp);
        }
      }
      vijk[((ZKPOP_N-1)*ZKPOP_TAU + j) * ZKPOP_M + k] = barrett_reduce(vk[k] - tmp);
    }
  }
  
  // compute commitments
  CALLOC8(hij, ZKPOP_N*ZKPOP_TAU*32);
  for (i = 0; i < ZKPOP_N; i++)
  {
    for (j = 0; j < ZKPOP_TAU; j++)
    {
      hash_h(&hij[(i*ZKPOP_TAU + j) * 32], (uint8_t*)&vijk[(i*ZKPOP_TAU + j) * ZKPOP_M], ZKPOP_M * sizeof(int16_t));
    }
  }
  
  // compute h
  hash_h(h, hij, ZKPOP_N*ZKPOP_TAU*32);
  
  // sample audit_bundles
  sample_audit_bundles(audit_bundles, h);
  
  // generate Sij, Eij
  CALLOCPOLYVEC(S, ZKPOP_N*ZKPOP_TAU);
  CALLOCPOLYVEC(E, ZKPOP_N*ZKPOP_TAU);
  for (i = 0; i < ZKPOP_N; i++)
  {
    for (j = 0; j < ZKPOP_TAU; j++)
    {
      count = 0;
      for (k = 0; k < ZKPOP_M; k++)
      {
        if (count < ZKPOP_M-ZKPOP_SIGMA && k == audit_bundles[count])
        {
          count += 1;
        } else {
          if (k - count < KYBER_K*KYBER_N)
          {
            S[i*ZKPOP_TAU + j].vec[(k-count)/KYBER_N].coeffs[(k-count)%KYBER_N] = vijk[(i*ZKPOP_TAU + j) * ZKPOP_M + k];
          } else {
            E[i*ZKPOP_TAU + j].vec[(k-count-KYBER_K*KYBER_N)/KYBER_N].coeffs[(k-count-KYBER_K*KYBER_N)%KYBER_N] = vijk[(i*ZKPOP_TAU + j) * ZKPOP_M + k];
          }
        }
      }
    }
  }
  
  // generate S_true and E_true
  count = 0;
  for (k = 0; k < ZKPOP_M; k++)
  {
    if (count < ZKPOP_M-ZKPOP_SIGMA && k == audit_bundles[count])
    {
      count += 1;
    } else {
      if (k - count < KYBER_K*KYBER_N)
      {
        S_true.vec[(k-count)/KYBER_N].coeffs[(k-count)%KYBER_N] = vk[k];
      } else {
        E_true.vec[(k-count-KYBER_K*KYBER_N)/KYBER_N].coeffs[(k-count-KYBER_K*KYBER_N)%KYBER_N] = vk[k];
      }
    }
  }
  
  // generate matrix A
  gen_a(a, publicseed);
  
  // compute B_true
  polyvec_ntt(&S_true);
  polyvec_ntt(&E_true);
  for(i=0; i < KYBER_K; i++) 
  {
    polyvec_basemul_acc_montgomery(&B_true.vec[i], &a[i], &S_true);
    poly_tomont(&B_true.vec[i]);
  }

  polyvec_add(&B_true, &B_true, &E_true);
  polyvec_reduce(&B_true);
  
  
  // compute Bij
  CALLOCPOLYVEC(B, ZKPOP_N*ZKPOP_TAU);
  CALLOC8(hBij, ZKPOP_N*ZKPOP_TAU*32); // potentially big
  // Bij = A*Sij + Eij, then commit to this value
  for (i = 0; i < ZKPOP_N; i++)
  {
    for (j = 0; j < ZKPOP_TAU; j++)
    {
      polyvec_ntt(&S[i*ZKPOP_TAU + j]);
      polyvec_ntt(&E[i*ZKPOP_TAU + j]);
      for(mi = 0; mi < KYBER_K; mi++) 
      {
        polyvec_basemul_acc_montgomery(&B[i*ZKPOP_TAU + j].vec[mi], &a[mi], &S[i*ZKPOP_TAU + j]);
        poly_tomont(&B[i*ZKPOP_TAU + j].vec[mi]);
      }

      polyvec_add(&B[i*ZKPOP_TAU + j], &B[i*ZKPOP_TAU + j], &E[i*ZKPOP_TAU + j]);
      polyvec_reduce(&B[i*ZKPOP_TAU + j]);
      
      // TODO do a polyvec_tobytes before to reduce the hash input size
      // this packed polyvec could be re-used later when packing the proof
      hash_h(&hBij[(i*ZKPOP_TAU + j)*32], (uint8_t*)&B[i*ZKPOP_TAU + j], sizeof(polyvec));
    }
  }
  
  // compute hB
  hash_h(hB, hBij, ZKPOP_N*ZKPOP_TAU*32);
  free(hBij);
  
  // sample audit_parties
  sample_audit_parties(audit_parties, h, hB, publicseed, &B_true);
  
  
  
  // pack everything into proof buffer
  *zkpop_size = proof_size = KYBER_ZKPOP_MAXBYTES;;
  CALLOC8(*zkpop, proof_size);
  
  uint8_t *zkpop_cur = *zkpop;
  STOREZKPOP(audit_parties, sizeof(uint16_t) * ZKPOP_TAU);
  STOREZKPOP(audit_bundles, sizeof(uint16_t) * (ZKPOP_M - ZKPOP_SIGMA));
  STOREZKPOP(h, 32);
  STOREZKPOP(hB, 32);
  
  // send seed tree nodes such that seedij[i][j] with i = audit_parties[j] can NOT be reconstructed
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (audit_parties[j] != ZKPOP_N-1)
    {
      // send all seeds but one
      i = ZKPOP_N + audit_parties[j] - 2;
      while (i != 0)
      {
        // add sibling of current node
        if (i & 1)
        {
          STOREZKPOP(&seedtree[(j*((ZKPOP_N-1)*2-1) + i+1)*32], 32);
        } else {
          STOREZKPOP(&seedtree[(j*((ZKPOP_N-1)*2-1) + i-1)*32], 32);
        }
        // go to parent node
        i = (i-1) / 2;
      }
    } else {
      STOREZKPOP(&seedtree[(j*((ZKPOP_N-1)*2-1))*32], 32);
    }
  }
  free(seedtree);
    
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (audit_parties[j] != ZKPOP_N-1)
    {
      for (k = 0; k < ZKPOP_M - ZKPOP_SIGMA; k++)
      {
        STOREZKPOP(&vijk[(audit_parties[j]*ZKPOP_TAU + j) * ZKPOP_M + audit_bundles[k]], sizeof(uint16_t));
      }
    }
  }
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    for (k = 0; k < ZKPOP_M - ZKPOP_SIGMA; k++)
    {
      STOREZKPOP(&vijk[((ZKPOP_N-1)*ZKPOP_TAU + j) * ZKPOP_M + audit_bundles[k]], sizeof(uint16_t));
    }
  }
  free(vijk);
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    STOREZKPOP(&hij[(ZKPOP_TAU*audit_parties[j] + j)*32], 32);
    if (audit_parties[j] != ZKPOP_N-1)
    {
      STOREZKPOP(&hij[(ZKPOP_TAU*(ZKPOP_N-1) + j)*32], 32);
    }
  }
  free(hij);
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    STOREZKPOP_POLYVEC(&B[audit_parties[j]*ZKPOP_TAU + j]);
  }
  free(B);
    
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (audit_parties[j] != ZKPOP_N-1)
    {
      STOREZKPOP_POLYVEC(&S[(ZKPOP_N-1)*ZKPOP_TAU + j]);
      STOREZKPOP_POLYVEC(&E[(ZKPOP_N-1)*ZKPOP_TAU + j]);
    }
  }
  free(S);
  free(E);
  
  *zkpop_size -= proof_size;
  
  pack_sk(sk, &S_true); // S
  pack_pk(pk, &B_true, publicseed); // B
  return 0;
}

static void build_subtree(uint8_t *seedtree, size_t start) // TODO iterative
{
  if (2*start+1 < (ZKPOP_N-1)*2-1)
  {
    prf(&seedtree[(2*start+1)*32], 64, &seedtree[start*32], 0);
    build_subtree(seedtree, 2*start+1);
    build_subtree(seedtree, 2*start+2);
  }
}

int crypto_zkpop_verify(const unsigned char *pk, const unsigned char *zkpop, unsigned long zkpop_size)
{
  size_t i,j,k,mi,count;
  uint16_t audit_parties[ZKPOP_TAU], audit_parties_check[ZKPOP_TAU];
  uint16_t audit_bundles[ZKPOP_M-ZKPOP_SIGMA], audit_bundles_check[ZKPOP_M-ZKPOP_SIGMA];
  const uint8_t *zkpop_cur = zkpop;
  int16_t vjk[ZKPOP_TAU][ZKPOP_M-ZKPOP_SIGMA] = {0};
  uint8_t *hij;
  int16_t *vijk;
  polyvec *S, *E, *B;
  uint8_t *hBij;
  polyvec B_true;
  uint8_t h[32], hB[32], h_check[32], hB_check[32], publicseed[KYBER_SYMBYTES];
  uint8_t *seedtree;
  polyvec a[KYBER_K];
  
  unpack_pk(&B_true, publicseed, pk);
  
  CALLOC16(vijk, ZKPOP_N*ZKPOP_TAU*ZKPOP_M);
  CALLOCPOLYVEC(E, ZKPOP_N*ZKPOP_TAU);
  CALLOCPOLYVEC(S, ZKPOP_N*ZKPOP_TAU);
  CALLOCPOLYVEC(B, ZKPOP_N*ZKPOP_TAU);
  
  // allocate memory that is potentially big in some parametrizations
  CALLOC8(seedtree, ZKPOP_TAU*((ZKPOP_N-1)*2-1)*32);
  CALLOC8(hij, ZKPOP_N*ZKPOP_TAU*32);
  
  LOADZKPOP(audit_parties, sizeof(uint16_t) * ZKPOP_TAU);
  for (j = 0; j < ZKPOP_TAU; j++) // check for valid range
  {
    if (audit_parties[j] >= ZKPOP_N)
    {
      return __LINE__;
    }
  }
  LOADZKPOP(audit_bundles, sizeof(uint16_t)*(ZKPOP_M-ZKPOP_SIGMA));
  for (k = 1; k < ZKPOP_M-ZKPOP_SIGMA; k++) // check for valid range and sortedness
  {
    if (audit_bundles[k-1] >= audit_bundles[k])
    {
      return __LINE__;
    }
    if (audit_bundles[k] >= ZKPOP_M)
    {
      return __LINE__;
    }
  }
  
  LOADZKPOP(h, 32);
  LOADZKPOP(hB, 32);
  
  // unpack seed tree nodes such that seedij[i][j] where i = audit_parties[j] can NOT be reconstructed
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (audit_parties[j] != ZKPOP_N-1)
    {
      i = ZKPOP_N + audit_parties[j] - 2;
      while (i != 0)
      {
        // load sibling of current node
        if (i & 1)
        {
          i += 1;
        } else {
          i -= 1;
        }
        LOADZKPOP(&seedtree[(j*((ZKPOP_N-1)*2-1) + i)*32], 32);
        // re-build subtree
        build_subtree(&seedtree[(j*((ZKPOP_N-1)*2-1) + 0)*32], i);
        // go to parent node
        i = (i-1) / 2;
      }
    } else {
      LOADZKPOP(&seedtree[(j*((ZKPOP_N-1)*2-1) + 0)*32], 32);
      build_subtree(&seedtree[(j*((ZKPOP_N-1)*2-1) + 0)*32], 0);
    }
    
    for (i = 0; i < ZKPOP_N-1; i++)
    {
      if (i != audit_parties[j])
      {
        // re-sample vijk
        uniform_vec(&vijk[(i*ZKPOP_TAU + j) * ZKPOP_M], ZKPOP_M, &seedtree[(j*((ZKPOP_N-1)*2-1) + ZKPOP_N + i - 2)*32]);
        // re-compute commitments
        hash_h(&hij[(i*ZKPOP_TAU + j)*32], (uint8_t*)(&vijk[(i*ZKPOP_TAU + j) * ZKPOP_M]), ZKPOP_M*sizeof(uint16_t));
      }
    }
  }
  free(seedtree);
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (audit_parties[j] != ZKPOP_N-1)
    {
      for (k = 0; k < ZKPOP_M - ZKPOP_SIGMA; k++)
      {
        LOADZKPOP(&vijk[(audit_parties[j]*ZKPOP_TAU + j) * ZKPOP_M + audit_bundles[k]], sizeof(uint16_t));
      }
    }
  }
  
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    for (k = 0; k < ZKPOP_M - ZKPOP_SIGMA; k++)
    {
      LOADZKPOP(&vijk[((ZKPOP_N-1)*ZKPOP_TAU + j) * ZKPOP_M + audit_bundles[k]], sizeof(uint16_t));
    }
  }
  
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    LOADZKPOP(&hij[(ZKPOP_TAU*audit_parties[j] + j)*32], 32);
    if (audit_parties[j] != ZKPOP_N-1)
    {
      LOADZKPOP(&hij[(ZKPOP_TAU*(ZKPOP_N-1) + j)*32], 32);
    }
  }
  
  CALLOC8(hBij, ZKPOP_N*ZKPOP_TAU*32); // potentially big
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    LOADZKPOP_POLYVEC(&B[audit_parties[j]*ZKPOP_TAU + j]);
    polyvec_reduce(&B[audit_parties[j]*ZKPOP_TAU + j]); // TODO conditional sub would do it as well
    // recompute hBij
    hash_h(&hBij[(audit_parties[j]*ZKPOP_TAU + j) * 32], (uint8_t*)&B[audit_parties[j]*ZKPOP_TAU + j], sizeof(polyvec));
  }
    
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    if (audit_parties[j] != ZKPOP_N-1)
    {
      LOADZKPOP_POLYVEC(&S[(ZKPOP_N-1)*ZKPOP_TAU + j]);
      LOADZKPOP_POLYVEC(&E[(ZKPOP_N-1)*ZKPOP_TAU + j]);
    }
  }
  if (zkpop_size != 0)
  {
    printf("remaining bytes: %lu\n", zkpop_size); // TODO remove this sanity check
    return __LINE__;
  }
  
  // generate Sij, Eij and re-compute audited vjk
  for (i = 0; i < ZKPOP_N-1; i++)
  {
    for (j = 0; j < ZKPOP_TAU; j++)
    {
      if (i != audit_parties[j])
      {
        count = 0;
        for (k = 0; k < ZKPOP_M; k++)
        {
          if (count < ZKPOP_M-ZKPOP_SIGMA && k == audit_bundles[count])
          {
            // TODO can we do lazy reduction here?
            vjk[j][count] = barrett_reduce(vjk[j][count] + vijk[(i*ZKPOP_TAU + j) * ZKPOP_M + k]);
            count += 1;
          } else {
            if (k-count < KYBER_K*KYBER_N)
            {
              S[i*ZKPOP_TAU + j].vec[(k-count)/KYBER_N].coeffs[(k-count)%KYBER_N] = vijk[(i*ZKPOP_TAU + j) * ZKPOP_M + k];
            } else {
              E[i*ZKPOP_TAU + j].vec[(k-count-KYBER_K*KYBER_N)/KYBER_N].coeffs[(k-count-KYBER_K*KYBER_N)%KYBER_N] = vijk[(i*ZKPOP_TAU + j) * ZKPOP_M + k];
            }
          }
        }
      } else {
        for (k = 0; k < ZKPOP_M - ZKPOP_SIGMA; k++)
        {
          // TODO lazy reduction?
          vjk[j][k] = barrett_reduce(vjk[j][k] + vijk[(i*ZKPOP_TAU + j) * ZKPOP_M + audit_bundles[k]]);
        }
      }
    }
  }
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    for (k = 0; k < ZKPOP_M-ZKPOP_SIGMA; k++)
    {
      vjk[j][k] = barrett_reduce(vjk[j][k] + vijk[((ZKPOP_N-1)*ZKPOP_TAU + j) * ZKPOP_M + audit_bundles[k]]);
    }
  }
  
  // check if all vjk are (1) small and (2) the same for each j
  for (k = 0; k < ZKPOP_M-ZKPOP_SIGMA; k++)
  {
    if (vjk[0][k] < -KYBER_ETA1 || vjk[0][k] > KYBER_ETA1)
    {
      return __LINE__;
    }
  }
  for (j = 1; j < ZKPOP_TAU; j++)
  {
    for (k = 0; k < ZKPOP_M-ZKPOP_SIGMA; k++)
    {
      if (vjk[j][k] != vjk[0][k])
      {
        return __LINE__;
      }
    }
  }
  // check h
  hash_h(h_check, hij, ZKPOP_TAU*ZKPOP_N*32);
  if (memcmp(h_check, h, 32) != 0)
  {
    return __LINE__;
  }
  
  // recompute audit_bundles
  sample_audit_bundles(audit_bundles_check, h_check);
  if (memcmp(audit_bundles, audit_bundles_check, sizeof(uint16_t)*(ZKPOP_M-ZKPOP_SIGMA)) != 0)
  {
    return __LINE__;
  }

  // generate matrix A
  gen_a(a, publicseed);
  
  // compute Bij for all i != r_j
  for (i = 0; i < ZKPOP_N; i++)
  {
    for (j = 0; j < ZKPOP_TAU; j++)
    {
      if (i != audit_parties[j])
      {
        if (i != ZKPOP_N-1)
        {
          polyvec_ntt(&S[i*ZKPOP_TAU + j]);
          polyvec_ntt(&E[i*ZKPOP_TAU + j]);
        }
        
        for(mi = 0; mi < KYBER_K; mi++) 
        {
          polyvec_basemul_acc_montgomery(&B[i*ZKPOP_TAU + j].vec[mi], &a[mi], &S[i*ZKPOP_TAU + j]);
          poly_tomont(&B[i*ZKPOP_TAU + j].vec[mi]);
        }

        polyvec_add(&B[i*ZKPOP_TAU + j], &B[i*ZKPOP_TAU + j], &E[i*ZKPOP_TAU + j]);
        polyvec_reduce(&B[i*ZKPOP_TAU + j]);
        
        // recompute commitment
        hash_h(&hBij[(i*ZKPOP_TAU + j)*32], (uint8_t*)&B[i*ZKPOP_TAU + j], sizeof(polyvec));
      }
    }
  }
  
  // check hB
  hash_h(hB_check, hBij, ZKPOP_N*ZKPOP_TAU*32);
  if (memcmp(hB, hB_check, 32) != 0)
  {
    return __LINE__;
  }
  free(hBij);
  
  // check if B = sum(Bij)
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    polyvec sum = {0};
    for (i = 0; i < ZKPOP_N; i++)
    {
      polyvec_add(&sum, &sum, &B[i*ZKPOP_TAU + j]);
      if ((i+1)%9 == 0)
      {
        polyvec_reduce(&sum);
      }
    }
    polyvec_reduce(&sum);
    polyvec_reduce(&B_true);
//     for (i = 0; i < KYBER_K; i++)
//     {
//       for (k = 0; k < KYBER_N; k++)
//       {
//         printf("%2lu, %3lu: %d, %d, %d\n", i,k, B_true.vec[i].coeffs[k], sum.vec[i].coeffs[k], B_true.vec[i].coeffs[k]-sum.vec[i].coeffs[k]);
//         
//       }
//     }
    if (memcmp(&sum, &B_true, sizeof(polyvec)) != 0)
    {
      return __LINE__;
    }
  }
  
  // recompute audit parties
  sample_audit_parties(audit_parties_check, h, hB, publicseed, &B_true);
  if (memcmp(audit_parties, audit_parties_check, sizeof(uint16_t)*ZKPOP_TAU) != 0)
  {
    return __LINE__;
  }
  
  free(hij);
  free(vijk);
  free(B);
  free(S);
  free(E);
  return 0;
}

/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
void indcpa_keypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  randombytes(buf, KYBER_SYMBYTES);
  hash_g(buf, buf, KYBER_SYMBYTES);

  gen_a(a, publicseed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&skpv.vec[i], noiseseed, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&e.vec[i], noiseseed, nonce++);

  polyvec_ntt(&skpv);
  polyvec_ntt(&e);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_basemul_acc_montgomery(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont(&pkpv.vec[i]);
  }

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c: pointer to output ciphertext
*                            (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m: pointer to input message
*                                  (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk: pointer to input public key
*                                   (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins used as seed
*                                      (of length KYBER_SYMBYTES) to deterministically
*                                      generate all randomness
**************************************************/
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], b;
  poly v, k, epp;

  unpack_pk(&pkpv, seed, pk);
  poly_frommsg(&k, m);
  gen_at(at, seed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(sp.vec+i, coins, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta2(ep.vec+i, coins, nonce++);
  poly_getnoise_eta2(&epp, coins, nonce++);

  polyvec_ntt(&sp);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_basemul_acc_montgomery(&b.vec[i], &at[i], &sp);

  polyvec_basemul_acc_montgomery(&v, &pkpv, &sp);

  polyvec_invntt_tomont(&b);
  poly_invntt_tomont(&v);

  polyvec_add(&b, &b, &ep);
  poly_add(&v, &v, &epp);
  poly_add(&v, &v, &k);
  polyvec_reduce(&b);
  poly_reduce(&v);

  pack_ciphertext(c, &b, &v);
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m: pointer to output decrypted message
*                            (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c: pointer to input ciphertext
*                                  (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec b, skpv;
  poly v, mp;

  unpack_ciphertext(&b, &v, c);
  unpack_sk(&skpv, sk);

  polyvec_ntt(&b);
  polyvec_basemul_acc_montgomery(&mp, &skpv, &b);
  poly_invntt_tomont(&mp);

  poly_sub(&mp, &v, &mp);
  poly_reduce(&mp);

  poly_tomsg(m, &mp);
}
