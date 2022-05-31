/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: functions for FrodoKEM-976
*           Instantiates "frodo_macrify.c" with the necessary matrix arithmetic functions
*********************************************************************************************/

#include "api_frodo976.h"
#include "frodo_macrify.h"


// Parameters for "FrodoKEM-976"
#define PARAMS_N 976
#define PARAMS_NBAR 8
#define PARAMS_LOGQ 16
#define PARAMS_Q (1 << PARAMS_LOGQ)
#define PARAMS_EXTRACTED_BITS 3
#define PARAMS_STRIPE_STEP 8
#define PARAMS_PARALLEL 4
#define PARAMS_BETA 10
#define BYTES_SEED_A 16
#define BYTES_MU (PARAMS_EXTRACTED_BITS*PARAMS_NBAR*PARAMS_NBAR)/8
#define BYTES_PKHASH CRYPTO_BYTES

#define ZKPOP_M 19485
#define ZKPOP_SIGMA 15616
#ifndef ZKPOP_N
#define ZKPOP_N 4
#endif
#ifndef ZKPOP_TAU
#define ZKPOP_TAU 96
#endif

#define ZKPOP_SYMBYTES 48

#if (PARAMS_NBAR % 8 != 0)
#error You have modified the cryptographic parameters. FrodoKEM assumes PARAMS_NBAR is a multiple of 8.
#endif

// Selecting SHAKE XOF function for the KEM and noise sampling
#define shake     shake256

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

// CDF table
uint16_t CDF_TABLE[11] = {5638, 15915, 23689, 28571, 31116, 32217, 32613, 32731, 32760, 32766, 32767};
uint16_t CDF_TABLE_LEN = 11;

#define crypto_kem_keypair            crypto_kem_keypair_Frodo976
#define crypto_kem_enc                crypto_kem_enc_Frodo976
#define crypto_kem_dec                crypto_kem_dec_Frodo976

#define crypto_kem_keypair_nizkpop    crypto_kem_keypair_nizkpop_Frodo976
#define crypto_nizkpop_verify         crypto_nizkpop_verify_Frodo976

#include "kem.c"
#include "noise.c"
#include "zkpop.c"
#if defined(USE_REFERENCE)
#include "frodo_macrify_reference.c"
#else
#include "frodo_macrify.c"
#endif
