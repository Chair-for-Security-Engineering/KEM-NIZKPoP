#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "log.h"

#if !defined(USE_AVX2)
    #include "sha3/fips202.h"
#else
    #include <immintrin.h>
    #include "sha3/fips202x4.h"
#endif

#if defined(USE_AES128_FOR_A)
    #include "aes/aes.h"
#endif    


#define NIZKPOP_MAX_BYTES (ZKPOP_SYMBYTES * (3 + ZKPOP_TAU * (LOG(ZKPOP_N) + 2)) + sizeof(uint16_t)*ZKPOP_M * ZKPOP_TAU + sizeof(uint16_t)*(ZKPOP_M-ZKPOP_SIGMA))

// TODO return codes instead of __LINE__


#define CALLOC16(x, size) uint16_t *x; do { x = calloc((size_t)size, sizeof(uint16_t)); if (x == NULL) return __LINE__; } while(0)
#define CALLOC8(x, size) uint8_t *x; do { x = calloc((size_t)size, sizeof(uint8_t)); if (x == NULL) return __LINE__; } while(0)
#define FREE(x) free(x)


#if ZKPOP_N > 65536
#error "Too many parties."
#elif ZKPOP_N > 256
typedef uint16_t party_t;
#else 
typedef uint8_t party_t;
#endif

// static void printhash(const char* str, uint8_t x[ZKPOP_SYMBYTES])
// {
//   printf("%s", str);
//   for (size_t i = 0; i < ZKPOP_SYMBYTES; i++)
//   {
//     printf("%02x", x[i]);
//   }
//   printf("\n");
// }

/**
 * Generate the complete matrix a.
 * 
 * As 
 */
static int frodo_expand_a(int16_t *a_row, const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) row-wise
    int i, k;
    
#if defined(USE_AES128_FOR_A)
    int16_t a_row_temp[4*PARAMS_N] = {0};                       // Take four lines of A at once       
#if defined(NO_OPENSSL)
    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);   
#else
    EVP_CIPHER_CTX *aes_key_schedule;    
    int len;
    if (!(aes_key_schedule = EVP_CIPHER_CTX_new())) handleErrors();    
    if (1 != EVP_EncryptInit_ex(aes_key_schedule, EVP_aes_128_ecb(), NULL, seed_A, NULL)) handleErrors();    
#endif
                                     
    for (int j = 0; j < PARAMS_N; j += PARAMS_STRIPE_STEP) {
        a_row_temp[j + 1 + 0*PARAMS_N] = UINT16_TO_LE(j);       // Loading values in the little-endian order
        a_row_temp[j + 1 + 1*PARAMS_N] = UINT16_TO_LE(j);
        a_row_temp[j + 1 + 2*PARAMS_N] = UINT16_TO_LE(j);
        a_row_temp[j + 1 + 3*PARAMS_N] = UINT16_TO_LE(j);
    }

    for (i = 0; i < PARAMS_N; i += 4) {
        for (int j = 0; j < PARAMS_N; j += PARAMS_STRIPE_STEP) {    // Go through A, four rows at a time
            a_row_temp[j + 0*PARAMS_N] = UINT16_TO_LE(i+0);     // Loading values in the little-endian order                                
            a_row_temp[j + 1*PARAMS_N] = UINT16_TO_LE(i+1);
            a_row_temp[j + 2*PARAMS_N] = UINT16_TO_LE(i+2);
            a_row_temp[j + 3*PARAMS_N] = UINT16_TO_LE(i+3);
        }

#if defined(NO_OPENSSL)
        AES128_ECB_enc_sch((uint8_t*)a_row_temp, 4*PARAMS_N*sizeof(int16_t), aes_key_schedule, (uint8_t*)(a_row + i*PARAMS_N));
#else   
        if (1 != EVP_EncryptUpdate(aes_key_schedule, (uint8_t*)a_row, &len, (uint8_t*)a_row_temp, 4*PARAMS_N*sizeof(int16_t))) handleErrors();
#endif
#elif defined (USE_SHAKE128_FOR_A)       
#if !defined(USE_AVX2)
    uint8_t seed_A_separated[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin = (uint16_t*)&seed_A_separated;
    memcpy(&seed_A_separated[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i += 4) {
        seed_A_origin[0] = UINT16_TO_LE(i + 0);
        shake128((unsigned char*)(a_row + (i+0)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 1);
        shake128((unsigned char*)(a_row + (i+1)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 2);
        shake128((unsigned char*)(a_row + (i+2)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 3);
        shake128((unsigned char*)(a_row + (i+3)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
#else
    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_2[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_3[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t* seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    uint16_t* seed_A_origin_2 = (uint16_t*)&seed_A_separated_2;
    uint16_t* seed_A_origin_3 = (uint16_t*)&seed_A_separated_3;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_2[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_3[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i += 4) {
        seed_A_origin_0[0] = UINT16_TO_LE(i + 0);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 1);
        seed_A_origin_2[0] = UINT16_TO_LE(i + 2);
        seed_A_origin_3[0] = UINT16_TO_LE(i + 3);
        shake128x4((unsigned char*)(a_row + (i+0)*PARAMS_N), (unsigned char*)(a_row + (i+1)*PARAMS_N), (unsigned char*)(a_row + (i+2)*PARAMS_N), (unsigned char*)(a_row + (i+3)*PARAMS_N), 
                    (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, seed_A_separated_2, seed_A_separated_3, 2 + BYTES_SEED_A);
#endif
#endif
        for (k = 0; k < 4 * PARAMS_N; k++) {
            a_row[k + i*PARAMS_N] = LE_TO_UINT16(a_row[k + i*PARAMS_N]);
        }
    }
    
#if defined(USE_AES128_FOR_A)
    AES128_free_schedule(aes_key_schedule);
#endif
    return 1;
}

static void build_subtree(uint8_t *seedtree, size_t start) // TODO iterative
{
  if (2*start+1 < ZKPOP_N*2-1)
  {
    SHAKE(&seedtree[(2*start+1)*ZKPOP_SYMBYTES], 2*ZKPOP_SYMBYTES, &seedtree[start*ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
    build_subtree(seedtree, 2*start+1);
    build_subtree(seedtree, 2*start+2);
  }
}

/*************************************************
* Name:        sample_audit_parties
*
* Description: Sample ZKPOP_TAU uniform random parties (range 0 to ZKPOP_N-1)
**************************************************/
static void sample_audit_parties (party_t audit_parties[ZKPOP_TAU], const uint8_t h2[ZKPOP_SYMBYTES])
{
  keccak_state state;
  SHAKE_absorb_once(&state, h2, ZKPOP_SYMBYTES);
/* ZKPOP_N might be a power of two (=no rejection sampling)
 * or not (=rejection sampling)
 */
#if ((ZKPOP_N-1)&ZKPOP_N) == 0 && ZKPOP_N <= 65536
  SHAKE_squeeze((uint8_t*)audit_parties, ZKPOP_TAU*sizeof(party_t), &state);
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
    for (i = 0; i < SHAKE_RATE/8; i++)
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
    for (i = 0; i < SHAKE_RATE/8; i += 2)
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

#define ZKPOP_M_MASK (ZKPOP_M | (ZKPOP_M>>1) | (ZKPOP_M>>2) | (ZKPOP_M>>3) | (ZKPOP_M>>4) | (ZKPOP_M>>5) | (ZKPOP_M>>6) | (ZKPOP_M>>7) | (ZKPOP_M>>8) | (ZKPOP_M>>9) | (ZKPOP_M>>10) | (ZKPOP_M>>11) | (ZKPOP_M>>12) | (ZKPOP_M>>13) | (ZKPOP_M>>14) | (ZKPOP_M>>15))
static void sample_audit_bundles (uint8_t audit_bundles[(ZKPOP_M+7)/8], const uint8_t h[ZKPOP_SYMBYTES])
{
#if ZKPOP_M > 65536
#error "sample_audit_bundles is implemented for ZKPOP_M <= 2^16"
#endif
  size_t count = ZKPOP_SIGMA,i,len=0;
  keccak_state state;
  uint8_t buf[SHAKE_RATE];
  SHAKE_absorb_once(&state, h, ZKPOP_SYMBYTES);
  do {
    SHAKE_squeezeblocks(buf, 1, &state);
    len = SHAKE_RATE;
    
    i = 0;
    // iterate over squeezed block (potentially drops some bits)
    while (len >= 2 && count < ZKPOP_M)
    {
      uint16_t tmp = (buf[i] | (buf[i+1] << 8)) & ZKPOP_M_MASK;
      if (tmp < count) // rejection sampling
      {
        audit_bundles[count/8] &= ~(1 << (count % 8));
        audit_bundles[count/8] |= ((audit_bundles[tmp / 8] >> (tmp % 8)) & 1) << (count % 8);
        audit_bundles[tmp / 8] |= 1 << (tmp % 8);
        
        count += 1;
      }
      len -= 2;
      i += 2;
    }
  } while(count < ZKPOP_M);
}

/*************************************************
* Name:        generate_seeds
*
* Description: Generate the seeds needed for generating the bundles vijk.
*
* Arguments:   - uint8_t *seedtree: pointer to the seed tree buffer
**************************************************/
static void generate_seeds(uint8_t *seedtree)
{
  size_t i, j, count;
#ifndef USE_AVX2
  uint8_t *parent, *children;
#else
  uint8_t *parent[4], *children[4];
#endif
  
  // expand initial seed to ZKPOP_TAU root seeds
  SHAKE(seedtree, ZKPOP_TAU*ZKPOP_SYMBYTES, seedtree, ZKPOP_SYMBYTES);
  
  // copy root seeds to correct position
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    memcpy(seedtree + (ZKPOP_TAU-j-1)*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES, seedtree + (ZKPOP_TAU-j-1)*ZKPOP_SYMBYTES, ZKPOP_SYMBYTES);
  }
  // build tree for each execution
#ifndef USE_AVX2
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    count = 1;
    uint8_t *parent, *children;
    memcpy(&seedtree[j*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES], rootseeds[j], ZKPOP_SYMBYTES);
    parent = &seedtree[j*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES];
    children = &seedtree[j*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES + ZKPOP_SYMBYTES];
    
    while (count < ZKPOP_N*2-1)
    {
      SHAKE(children, 2*ZKPOP_SYMBYTES, parent, ZKPOP_SYMBYTES);
      children += 2*ZKPOP_SYMBYTES;
      parent += ZKPOP_SYMBYTES;
      count += 2;
    }
  }
#else
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
      SHAKEx4(children[0],
                 children[1],
                 children[2],
                 children[3],
                 2 * ZKPOP_SYMBYTES,
                 parent[0],
                 parent[1],
                 parent[2],
                 parent[3],
                 ZKPOP_SYMBYTES);
      for (i = 0; i < 4; i++)
      {
        children[i] += 2 * ZKPOP_SYMBYTES;
        parent[i] += ZKPOP_SYMBYTES;
      }
      count += 2;
    }
  }
#endif
}

static void generate_S_E(uint16_t S[2*PARAMS_N*PARAMS_NBAR], const uint16_t vk[ZKPOP_M], const uint8_t audit_bundles[(ZKPOP_M+7)/8])
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
      S[k-count] = vk[k];
    }
    audit_flags >>= 1;
  }
}

static void generate_S_E_openbundle(uint16_t S[2*PARAMS_N*PARAMS_NBAR], uint16_t openbundle[ZKPOP_M-ZKPOP_SIGMA], const uint16_t vk[ZKPOP_M], const uint8_t audit_bundles[(ZKPOP_M+7)/8])
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
      openbundle[count] = vk[k] % PARAMS_Q;
      count += 1;
    } else {
      S[k-count] = vk[k] % PARAMS_Q;
    }
    audit_flags >>= 1;
  }
}

#define STOREZKPOP(src, size) do {if (proof_size < size) { return 0; } memcpy(zkpop_cur, src, size); zkpop_cur += size; proof_size -= size;} while(0)
/*************************************************
* Name:        pack_nizkpop
*
* Description: Pack all values to the NIZKPoP buffer.
**************************************************/
static size_t pack_nizkpop(uint8_t *zkpop, 
                           const uint8_t h1[ZKPOP_SYMBYTES], 
                           const uint8_t h2[ZKPOP_SYMBYTES], 
                           const uint8_t salt[ZKPOP_SYMBYTES], 
                           const uint8_t seedtree[(ZKPOP_N*2-1)*ZKPOP_TAU*ZKPOP_SYMBYTES], 
                           const uint8_t *comij,
                           const uint16_t vk[ZKPOP_M],
                           const party_t audit_parties[ZKPOP_TAU], 
                           const uint8_t audit_bundles[(ZKPOP_M+7)/8])
{
  size_t i, j, k;
  uint8_t audit_flags = 0;
  uint8_t *zkpop_cur = zkpop;
  size_t proof_size = NIZKPOP_MAX_BYTES;
  
  STOREZKPOP(h1, ZKPOP_SYMBYTES);
  STOREZKPOP(h2, ZKPOP_SYMBYTES);
  STOREZKPOP(salt, ZKPOP_SYMBYTES);
  
  // store comij unaudited
  for (j = 0; j < ZKPOP_TAU; j++)
  {
    STOREZKPOP(&comij[(audit_parties[j]*ZKPOP_TAU + j) * ZKPOP_SYMBYTES], ZKPOP_SYMBYTES);
  }
  
  // store all deltabjk
  //STOREZKPOP(deltabjk[0], ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t));
  // this is already done
  if (proof_size >= ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t))
  {
    zkpop_cur += ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t);
    proof_size -= ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t);
  } else {
    return __LINE__;
  }
    
  // store all unaudited vk
  for (k = 0; k < ZKPOP_M; k++)
  {
    if ((k%8) == 0)
    {
      audit_flags = audit_bundles[k/8];
    }
    if (audit_flags&1)
    {
      STOREZKPOP(&vk[k], sizeof(uint16_t));
    } 
    audit_flags >>= 1;
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
  
  return NIZKPOP_MAX_BYTES - proof_size;
}

#define LOADZKPOP(dest, size) do {if (zkpop_size < size) { /*printf("bad zkpop: size\n");*/ return -1;} memcpy(dest, zkpop_cur, size); zkpop_cur += size; zkpop_size -= size;} while(0)
/*************************************************
* Name:        unpack_nizkpop
*
* Description: Unpack given NIZKPoP buffer.
* 
* Returns 0 if ok.
**************************************************/
static int unpack_nizkpop(party_t audit_parties[ZKPOP_TAU],
                          uint8_t audit_bundles[(ZKPOP_M+7)/8],
                          uint16_t vk[ZKPOP_M-ZKPOP_SIGMA],
                          uint8_t seedtree[(ZKPOP_N*2-1)*ZKPOP_TAU*ZKPOP_SYMBYTES],
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
  // this is done by assigning a pointer
  if (zkpop_size >= ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t))
  {
    zkpop_cur += ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t);
    zkpop_size -= ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t);
  } else {
    return __LINE__;
  }
  
  // unpack all audited vk
  LOADZKPOP(vk, sizeof(uint16_t)*(ZKPOP_M-ZKPOP_SIGMA));
  
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



int crypto_kem_keypair_nizkpop(unsigned char* pk, unsigned char* sk, unsigned char **zkpop, unsigned long *zkpop_size)
{
  // iterators
  size_t i, j, k, mi, ij;
  
  // ZKPOP values
  uint8_t h1[ZKPOP_SYMBYTES], h2[ZKPOP_SYMBYTES];
  uint8_t audit_bundles[(ZKPOP_M+7)/8] = {0};
  party_t audit_parties[ZKPOP_TAU];
  CALLOC8(seedtree, ZKPOP_TAU*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES);
#ifndef NO_RESAMPLING
  uint16_t bijk[4][ZKPOP_M];
#else
  CALLOC16(bijk, ZKPOP_N*ZKPOP_TAU*ZKPOP_M);
#endif
  uint16_t vk[ZKPOP_M], *deltabjk;
  CALLOC8(comij, ZKPOP_TAU*ZKPOP_N*ZKPOP_SYMBYTES);
  
  // ZKPoP buffers
  uint8_t comijbuf[4][ZKPOP_SYMBYTES*2 + 5];
  uint8_t pvhash[4][ZKPOP_SYMBYTES];
  uint16_t openbundletmp[4][ZKPOP_M-ZKPOP_SIGMA];
  keccak_state state;
  
  // Frodo values
  uint8_t randomness[ZKPOP_SYMBYTES + ZKPOP_SYMBYTES + BYTES_SEED_A];
  uint8_t *salt = randomness + ZKPOP_SYMBYTES; // ZKPoP salt
  uint8_t *seedrand = randomness + ZKPOP_SYMBYTES + ZKPOP_SYMBYTES;
  uint16_t S[2*PARAMS_N*PARAMS_NBAR] = {0};
  uint16_t B[4][PARAMS_N*PARAMS_NBAR] = {0};
  uint8_t *pk_seedA = &pk[0];
  ALIGN_HEADER(32) int16_t a_row[PARAMS_N*PARAMS_N] ALIGN_FOOTER(32) = {0};
  
  *zkpop_size = NIZKPOP_MAX_BYTES;
  *zkpop = calloc(NIZKPOP_MAX_BYTES, 1);
  if (*zkpop == NULL)
  {
    return __LINE__;
  }
  deltabjk = (uint16_t*)(*zkpop + (3+ZKPOP_TAU)*ZKPOP_SYMBYTES);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Step 1
  randombytes(randomness, ZKPOP_SYMBYTES + ZKPOP_SYMBYTES + BYTES_SEED_A);
  // generate vk
  SHAKE((uint8_t*)vk, ZKPOP_M*2, randomness, ZKPOP_SYMBYTES);
  frodo_sample_n(vk, ZKPOP_M);
  
  // generate seed tree
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
#ifndef NO_RESAMPLING
    SHAKEx4( 
      (uint8_t*)bijk[0], 
      (uint8_t*)bijk[1], 
      (uint8_t*)bijk[2], 
      (uint8_t*)bijk[3], 
      ZKPOP_M*sizeof(uint16_t),
      comijbuf[0],
      comijbuf[1],
      comijbuf[2],
      comijbuf[3],
      ZKPOP_SYMBYTES*2 + 5
    );
#else
    SHAKEx4( 
      (uint8_t*)&bijk[(ZKPOP_N * J[0] + I[0]) * ZKPOP_M], 
      (uint8_t*)&bijk[(ZKPOP_N * J[1] + I[1]) * ZKPOP_M], 
      (uint8_t*)&bijk[(ZKPOP_N * J[2] + I[2]) * ZKPOP_M], 
      (uint8_t*)&bijk[(ZKPOP_N * J[3] + I[3]) * ZKPOP_M], 
      ZKPOP_M*sizeof(uint16_t),
      comijbuf[0],
      comijbuf[1],
      comijbuf[2],
      comijbuf[3],
      ZKPOP_SYMBYTES*2 + 5
    );
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
      ZKPOP_SYMBYTES*2+5
    );
    
    // update deltabjk
    for (mi = 0; mi < nj; mi++)
    {
      for (k = 0; k < ZKPOP_M; k++)
      {
#ifndef NO_RESAMPLING
        deltabjk[J[mi] * ZKPOP_M + k] += bijk[mi][k];
#else
        deltabjk[J[mi] * ZKPOP_M + k] += bijk[(ZKPOP_N * J[mi] + I[mi]) * ZKPOP_M + k];
#endif
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
    for (k = 0; k < ZKPOP_M; k++)
    {
      deltabjk[j * ZKPOP_M + k] = (vk[k] - deltabjk[j * ZKPOP_M + k]) % PARAMS_Q;
    }
    SHAKE_absorb(&state, (uint8_t*)&deltabjk[j * ZKPOP_M], ZKPOP_M*sizeof(uint16_t));
  }
  
  // TODO add absorption of attrs
  SHAKE_finalize(&state);
  SHAKE_squeeze(h1, ZKPOP_SYMBYTES, &state);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Step 2
  // let V choose ZKPOP_SIGMA bundles from knowing h
  sample_audit_bundles(audit_bundles, h1);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Step 3
  // generate pk_seedA
  SHAKE(pk_seedA, BYTES_SEED_A, seedrand, BYTES_SEED_A);
  
  // generate A
  frodo_expand_a(a_row, pk_seedA);
  
  
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
#ifndef NO_RESAMPLING
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
    
    SHAKEx4( 
      (uint8_t*)bijk[0], 
      (uint8_t*)bijk[1], 
      (uint8_t*)bijk[2], 
      (uint8_t*)bijk[3], 
      ZKPOP_M*sizeof(uint16_t),
      comijbuf[0],
      comijbuf[1],
      comijbuf[2],
      comijbuf[3],
      ZKPOP_SYMBYTES*2 + 5
    );
#else
    for (mi = 0; mi < ijcnt; mi++)
    {
      I[mi] = (ij+mi)%ZKPOP_N;
      J[mi] = (ij+mi)/ZKPOP_N;
    }
#endif
    
    for (mi = 0; mi < ijcnt; mi++)
    {
      // update first party
      if (I[mi] == 0)
      {
        for (k = 0; k < ZKPOP_M; k++)
        {
#ifndef NO_RESAMPLING
          bijk[mi][k] += deltabjk[J[mi] * ZKPOP_M + k];
          bijk[mi][k] %= PARAMS_Q;
#else 
          bijk[(ZKPOP_N * J[mi] + I[mi]) * ZKPOP_M + k] += deltabjk[J[mi] * ZKPOP_M + k];
          bijk[(ZKPOP_N * J[mi] + I[mi]) * ZKPOP_M + k] %= PARAMS_Q;
#endif
        }
      }
      
      // generate Sij and Eij from bundle
#ifndef NO_RESAMPLING
      generate_S_E_openbundle(S, openbundletmp[mi], bijk[mi], audit_bundles);
#else
      generate_S_E_openbundle(S, openbundletmp[mi], &bijk[(ZKPOP_N * J[mi] + I[mi]) * ZKPOP_M], audit_bundles);
#endif
      
      for (i = 0; i < (PARAMS_N*PARAMS_NBAR); i++) 
      {    
          B[mi][i] = S[PARAMS_N*PARAMS_NBAR + i];
      }
      for (i = 0; i < PARAMS_N; i++)
      {
        for (k = 0; k < PARAMS_NBAR; k++)
        {
          uint16_t sum = 0;
          for (j = 0; j < PARAMS_N; j++)
          {
            sum += a_row[PARAMS_N*i + j] * S[k*PARAMS_N + j];
          }
          B[mi][PARAMS_NBAR*i + k] += sum;
          B[mi][PARAMS_NBAR*i + k] %= PARAMS_Q;
        }
      }
    }
    SHAKEx4(
      pvhash[0], 
      pvhash[1],
      pvhash[2],
      pvhash[3],
      ZKPOP_SYMBYTES,
      (uint8_t*)B[0],
      (uint8_t*)B[1],
      (uint8_t*)B[2],
      (uint8_t*)B[3],
      PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    SHAKE_absorb(&state, pvhash[0], ZKPOP_SYMBYTES*ijcnt);
    
    SHAKEx4(
      pvhash[0], 
      pvhash[1],
      pvhash[2],
      pvhash[3],
      ZKPOP_SYMBYTES,
      (uint8_t*)openbundletmp[0],
      (uint8_t*)openbundletmp[1],
      (uint8_t*)openbundletmp[2],
      (uint8_t*)openbundletmp[3],
      (ZKPOP_M-ZKPOP_SIGMA)*sizeof(uint16_t));
    SHAKE_absorb(&bundlestate, pvhash[0], ZKPOP_SYMBYTES*ijcnt);
  }
  
  // generate S and E from those vk values that are NOT selected for auditing
  generate_S_E(S, vk, audit_bundles);
  
  // compute B as in original key gen
  for (i = 0; i < (PARAMS_N*PARAMS_NBAR); i++) 
  {    
      B[0][i] = S[PARAMS_N*PARAMS_NBAR + i];
  }
  for (i = 0; i < PARAMS_N; i++)
  {
    for (k = 0; k < PARAMS_NBAR; k++)
    {
      uint16_t sum = 0;
      for (j = 0; j < PARAMS_N; j++)
      {
        sum += a_row[PARAMS_N*i + j] * S[k*PARAMS_N + j];
      }
      B[0][PARAMS_NBAR*i + k] += sum;
      B[0][PARAMS_NBAR*i + k] %= PARAMS_Q;
    }
  }
  
  // hash result
  SHAKE_absorb(&state, (uint8_t*)B[0], PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
  SHAKE_finalize(&bundlestate);
  SHAKE_squeeze(h2, ZKPOP_SYMBYTES, &bundlestate); // use h2 to finalize the audited bundles hash
  SHAKE_absorb(&state, h2, ZKPOP_SYMBYTES);
  SHAKE_finalize(&state);
  SHAKE_squeeze(h2, ZKPOP_SYMBYTES, &state);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Step 4

  // sample one hidden party per j in {1,...,tau}
  sample_audit_parties(audit_parties, h2);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Step 5
  // pack everything into proof buffer:
  *zkpop_size = pack_nizkpop(
    *zkpop, 
    h1, 
    h2,
    salt,
    seedtree, 
    comij,
    vk, 
    audit_parties, 
    audit_bundles
  );
  if (*zkpop_size == 0)
  {
    return __LINE__;
  }
#ifdef NO_RESAMPLING
  FREE(bijk);
#endif
  FREE(seedtree);
  FREE(comij);

  // do what is left from FrodoKEM's key generation
  // Outputs: public key pk (               BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 bytes)
  //          secret key sk (CRYPTO_BYTES + BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + 2*PARAMS_N*PARAMS_NBAR + BYTES_PKHASH bytes)
  uint8_t *pk_b = &pk[BYTES_SEED_A];
  uint8_t *sk_s = &sk[0];
  uint8_t *sk_pk = &sk[CRYPTO_BYTES];
  uint8_t *sk_S = &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES];
  uint8_t *sk_pkh = &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES + 2*PARAMS_N*PARAMS_NBAR];

  // Encode the second part of the public key
  frodo_pack(pk_b, CRYPTO_PUBLICKEYBYTES - BYTES_SEED_A, B[0], PARAMS_N*PARAMS_NBAR, PARAMS_LOGQ);

  // Add s, pk and S to the secret key
  randombytes(sk_s, CRYPTO_BYTES);
  memcpy(sk_pk, pk, CRYPTO_PUBLICKEYBYTES);
  for (size_t i = 0; i < PARAMS_N * PARAMS_NBAR; i++) {
      S[i] = UINT16_TO_LE(S[i]);
  }
  memcpy(sk_S, S, 2*PARAMS_N*PARAMS_NBAR);

  // Add H(pk) to the secret key
  shake(sk_pkh, BYTES_PKHASH, pk, CRYPTO_PUBLICKEYBYTES);
  return 0;
}

int crypto_nizkpop_verify(const unsigned char *pk, const unsigned char *zkpop, unsigned long zkpop_size)
{
  // iterators
  size_t i,j,k,mi,mul_i,mul_j,mul_k,i2;
  size_t count = 0, I[4] = {0};
  
  // ZKPoP values and buffers
  party_t audit_parties[ZKPOP_TAU] = {0};
  uint8_t audit_bundles[(ZKPOP_M+7)/8] = {0};
  uint8_t h1[ZKPOP_SYMBYTES], h2[ZKPOP_SYMBYTES], salt[ZKPOP_SYMBYTES], h1_check[ZKPOP_SYMBYTES], h2_check[ZKPOP_SYMBYTES];
  CALLOC8(seedtree, ZKPOP_TAU*(ZKPOP_N*2-1)*ZKPOP_SYMBYTES);
  uint16_t vk[ZKPOP_M-ZKPOP_SIGMA], missingshare[ZKPOP_M-ZKPOP_SIGMA];
  uint16_t bijk[4][ZKPOP_M] = {0}, *deltabjk;
  uint16_t Sij[2*PARAMS_N*PARAMS_NBAR];
  uint16_t Bij[4][PARAMS_N*PARAMS_NBAR], missingBij[PARAMS_N*PARAMS_NBAR];
  uint8_t Bijhash[ZKPOP_N][ZKPOP_SYMBYTES], openbundlehash[ZKPOP_N][ZKPOP_SYMBYTES], tmphash[4][ZKPOP_SYMBYTES];
  uint8_t comijbuf[4][ZKPOP_SYMBYTES*2 + 5] = {0};
  uint16_t openbundlebuf[4][ZKPOP_M-ZKPOP_SIGMA];
  keccak_state h1state, h2state, bundlestate;
  const uint8_t *zkpptr = zkpop + 3*ZKPOP_SYMBYTES; // used for absorbing the transmitted comij
  
  // Frodo values
  ALIGN_HEADER(32) int16_t a_row[PARAMS_N*PARAMS_N] ALIGN_FOOTER(32) = {0};
  uint16_t B[PARAMS_N*PARAMS_NBAR];
  
  int rv;
    
  if ((rv = unpack_nizkpop(
    audit_parties, 
    audit_bundles, 
    vk, 
    seedtree,
    salt,
    h2,
    h1,
    zkpop,
    zkpop_size
  )) != 0)
  {
    FREE(seedtree);
    return rv;
  }
  
  deltabjk = (uint16_t*)(zkpop + ZKPOP_SYMBYTES*(ZKPOP_TAU + 3));
  
  // generate A
  frodo_expand_a(a_row, pk);
  
  // unpack B
  frodo_unpack(B, PARAMS_N*PARAMS_NBAR, &pk[BYTES_SEED_A], CRYPTO_PUBLICKEYBYTES - BYTES_SEED_A, PARAMS_LOGQ);
  
  

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
    memcpy(missingBij, B, sizeof(uint16_t)*PARAMS_N*PARAMS_NBAR);
    memcpy(missingshare, vk, sizeof(uint16_t)*(ZKPOP_M-ZKPOP_SIGMA));
    
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
        SHAKEx4( 
          (uint8_t*)bijk[0], 
          (uint8_t*)bijk[1], 
          (uint8_t*)bijk[2], 
          (uint8_t*)bijk[3], 
          ZKPOP_M*sizeof(uint16_t),
          comijbuf[0],
          comijbuf[1],
          comijbuf[2],
          comijbuf[3],
          ZKPOP_SYMBYTES*2 + 5
        );
        
        for (mi = 0; mi < count; mi++)
        {
          // update first party
          if (I[mi] == 0)
          {
            for (k = 0; k < ZKPOP_M; k++)
            {
              bijk[mi][k] += deltabjk[j * ZKPOP_M + k];
              bijk[mi][k] %= PARAMS_Q;
            }
          }
          
          // generate Sij and Eij from bundle
          generate_S_E_openbundle(Sij, openbundlebuf[mi], bijk[mi], audit_bundles);
          
          for (mul_i = 0; mul_i < (PARAMS_N*PARAMS_NBAR); mul_i++) 
          {    
              Bij[mi][mul_i] = Sij[PARAMS_N*PARAMS_NBAR + mul_i];
          }
          for (mul_i = 0; mul_i < PARAMS_N; mul_i++)
          {
            for (mul_k = 0; mul_k < PARAMS_NBAR; mul_k++)
            {
              uint16_t sum = 0;
              for (mul_j = 0; mul_j < PARAMS_N; mul_j++)
              {
                sum += a_row[PARAMS_N*mul_i + mul_j] * Sij[mul_k*PARAMS_N + mul_j];
              }
              Bij[mi][PARAMS_NBAR*mul_i + mul_k] += sum;
              Bij[mi][PARAMS_NBAR*mul_i + mul_k] %= PARAMS_Q;
            }
          }
          
          for (k = 0; k < PARAMS_N*PARAMS_NBAR; k++)
          {
            missingBij[k] -= Bij[mi][k];
            missingBij[k] %= PARAMS_Q;
          }
          
          // update missingshare
          for (k = 0; k < ZKPOP_M-ZKPOP_SIGMA; k++)
          {
            missingshare[k] -= openbundlebuf[mi][k];
            missingshare[k] %= PARAMS_Q;
          }
        }
        
        SHAKEx4(
          tmphash[0],
          tmphash[1],
          tmphash[2],
          tmphash[3],
          ZKPOP_SYMBYTES,
          (uint8_t*)Bij[0],
          (uint8_t*)Bij[1],
          (uint8_t*)Bij[2],
          (uint8_t*)Bij[3],
          sizeof(uint16_t)*PARAMS_N*PARAMS_NBAR
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
          (uint8_t*)openbundlebuf[0],
          (uint8_t*)openbundlebuf[1],
          (uint8_t*)openbundlebuf[2],
          (uint8_t*)openbundlebuf[3],
          sizeof(uint16_t)*(ZKPOP_M-ZKPOP_SIGMA)
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
    for (k = 0; k < PARAMS_N*PARAMS_NBAR; k++)
    {
      missingBij[k] %= PARAMS_Q;
    }
    SHAKE(Bijhash[audit_parties[j]], ZKPOP_SYMBYTES, (uint8_t*)missingBij, sizeof(uint16_t)*PARAMS_N*PARAMS_NBAR);
    SHAKE_absorb(&h2state, Bijhash[0], ZKPOP_SYMBYTES*ZKPOP_N);
    
    for (k = 0; k < ZKPOP_M-ZKPOP_SIGMA; k++)
    {
      missingshare[k] %= PARAMS_Q;
    }
    SHAKE(openbundlehash[audit_parties[j]], ZKPOP_SYMBYTES, (uint8_t*)missingshare, sizeof(uint16_t)*(ZKPOP_M-ZKPOP_SIGMA));
    SHAKE_absorb(&bundlestate, openbundlehash[0], ZKPOP_SYMBYTES*ZKPOP_N);
  }
  SHAKE_absorb(&h2state, (uint8_t*)B, sizeof(uint16_t)*PARAMS_N*PARAMS_NBAR);
  
  SHAKE_finalize(&bundlestate);
  SHAKE_squeeze(tmphash[0], ZKPOP_SYMBYTES, &bundlestate);
  SHAKE_absorb(&h2state, tmphash[0], ZKPOP_SYMBYTES);
  SHAKE_finalize(&h2state);
  SHAKE_squeeze(h2_check, ZKPOP_SYMBYTES, &h2state);
  
  if (memcmp(h2, h2_check, ZKPOP_SYMBYTES) != 0)
  {
    FREE(seedtree);
    return __LINE__;
  }
  
  // finalize h1
  SHAKE_absorb(&h1state, (uint8_t*)&deltabjk[0], ZKPOP_TAU*ZKPOP_M*sizeof(uint16_t));
  SHAKE_finalize(&h1state);
  SHAKE_squeeze(h1_check, ZKPOP_SYMBYTES, &h1state);
  
  if (memcmp(h1, h1_check, ZKPOP_SYMBYTES) != 0)
  {
    FREE(seedtree);
    return __LINE__;
  }
  FREE(seedtree);
  return 0;
}
