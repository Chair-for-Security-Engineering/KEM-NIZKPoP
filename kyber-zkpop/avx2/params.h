#ifndef PARAMS_H
#define PARAMS_H

#include "log.h"

#ifndef KYBER_K
#define KYBER_K 3	/* Change this for different security strengths */
#endif

//#define KYBER_90S	/* Uncomment this if you want the 90S variant */

/* Don't change parameters below this line */
#if   (KYBER_K == 2)
#ifdef KYBER_90S
#define KYBER_NAMESPACE(s) pqcrystals_kyber512_90s_avx2_##s
#else
#define KYBER_NAMESPACE(s) pqcrystals_kyber512_avx2_##s
#endif
#elif (KYBER_K == 3)
#ifdef KYBER_90S
#define KYBER_NAMESPACE(s) pqcrystals_kyber768_90s_avx2_##s
#else
#define KYBER_NAMESPACE(s) pqcrystals_kyber768_avx2_##s
#endif
#elif (KYBER_K == 4)
#ifdef KYBER_90S
#define KYBER_NAMESPACE(s) pqcrystals_kyber1024_90s_avx2_##s
#else
#define KYBER_NAMESPACE(s) pqcrystals_kyber1024_avx2_##s
#endif
#else
#error "KYBER_K must be in {2,3,4}"
#endif

#define KYBER_N 256
#define KYBER_Q 3329

#define KYBER_SYMBYTES 32   /* size in bytes of hashes, and seeds */
#define KYBER_SSBYTES  32   /* size in bytes of shared key */

#define KYBER_POLYBYTES		384
#define KYBER_POLYVECBYTES	(KYBER_K * KYBER_POLYBYTES)

#if KYBER_K == 2
#define KYBER_ETA1 3
#define KYBER_POLYCOMPRESSEDBYTES    128
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 320)
#define ZKPOP_M 1280
#define ZKPOP_SYMBYTES 32
#ifndef ZKPOP_N
#define ZKPOP_N 4
#endif
#ifndef ZKPOP_TAU
#define ZKPOP_TAU 64
#endif
#elif KYBER_K == 3
#define KYBER_ETA1 2
#define KYBER_POLYCOMPRESSEDBYTES    128
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 320)
#define ZKPOP_M 1870
#define ZKPOP_SYMBYTES 48
#ifndef ZKPOP_N
#define ZKPOP_N 4
#endif
#ifndef ZKPOP_TAU
#define ZKPOP_TAU 96
#endif
#elif KYBER_K == 4
#define KYBER_ETA1 2
#define KYBER_POLYCOMPRESSEDBYTES    160
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 352)
#define ZKPOP_M 2493
#define ZKPOP_SYMBYTES 64
#ifndef ZKPOP_N
#define ZKPOP_N 4
#endif
#ifndef ZKPOP_TAU
#define ZKPOP_TAU 128
#endif
#endif

#define ZKPOP_SIGMA (KYBER_K * KYBER_N * 2)
#define ZKPOP_M_PACKEDBYTES ((ZKPOP_M/KYBER_N)*KYBER_POLYBYTES + sizeof(int16_t)*(ZKPOP_M%KYBER_N))
#define ZKPOP_M_SIGMA_PACKEDBYTES (((ZKPOP_M-ZKPOP_SIGMA)/KYBER_N)*KYBER_POLYBYTES + sizeof(int16_t)*((ZKPOP_M-ZKPOP_SIGMA)%KYBER_N))

#define KYBER_ZKPOP_MAXBYTES (ZKPOP_SYMBYTES * (3 + ZKPOP_TAU * (LOG(ZKPOP_N) + 2)) + ZKPOP_M_PACKEDBYTES * ZKPOP_TAU + ZKPOP_M_SIGMA_PACKEDBYTES)

#define KYBER_ETA2 2

#define KYBER_INDCPA_MSGBYTES       (KYBER_SYMBYTES)
#define KYBER_INDCPA_PUBLICKEYBYTES (KYBER_POLYVECBYTES + KYBER_SYMBYTES)
#define KYBER_INDCPA_SECRETKEYBYTES (KYBER_POLYVECBYTES)
#define KYBER_INDCPA_BYTES          (KYBER_POLYVECCOMPRESSEDBYTES + KYBER_POLYCOMPRESSEDBYTES)

#define KYBER_PUBLICKEYBYTES  (KYBER_INDCPA_PUBLICKEYBYTES)
/* 32 bytes of additional space to save H(pk) */
#define KYBER_SECRETKEYBYTES  (KYBER_INDCPA_SECRETKEYBYTES + KYBER_INDCPA_PUBLICKEYBYTES + 2*KYBER_SYMBYTES)
#define KYBER_CIPHERTEXTBYTES (KYBER_INDCPA_BYTES)

#endif
