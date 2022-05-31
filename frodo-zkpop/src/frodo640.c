/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: functions for FrodoKEM-640
*           Instantiates "frodo_macrify.c" with the necessary matrix arithmetic functions
*********************************************************************************************/

#include "api_frodo640.h"
#include "frodo_macrify.h"


// Parameters for "FrodoKEM-640"
#define PARAMS_N 640
#define PARAMS_NBAR 8
#define PARAMS_LOGQ 15
#define PARAMS_Q (1 << PARAMS_LOGQ)
#define PARAMS_EXTRACTED_BITS 2
#define PARAMS_STRIPE_STEP 8
#define PARAMS_PARALLEL 4
#define PARAMS_BETA 12
#define BYTES_SEED_A 16
#define BYTES_MU (PARAMS_EXTRACTED_BITS*PARAMS_NBAR*PARAMS_NBAR)/8
#define BYTES_PKHASH CRYPTO_BYTES

// ZKPoP parameters for Frodo 640
#define ZKPOP_M 13233
#define ZKPOP_SIGMA 10240
#ifndef ZKPOP_N
#define ZKPOP_N 4
#endif
#ifndef ZKPOP_TAU
#define ZKPOP_TAU 64
#endif

#define ZKPOP_SYMBYTES 32

#if (PARAMS_NBAR % 8 != 0)
#error You have modified the cryptographic parameters. FrodoKEM assumes PARAMS_NBAR is a multiple of 8.
#endif

// Selecting SHAKE XOF function for the KEM and noise sampling
#define shake     shake128

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

// CDF table
uint16_t CDF_TABLE[13] = {4643, 13363, 20579, 25843, 29227, 31145, 32103, 32525, 32689, 32745, 32762, 32766, 32767};
uint16_t CDF_TABLE_LEN = 13;

#define crypto_kem_keypair            crypto_kem_keypair_Frodo640
#define crypto_kem_enc                crypto_kem_enc_Frodo640
#define crypto_kem_dec                crypto_kem_dec_Frodo640

#define crypto_kem_keypair_nizkpop    crypto_kem_keypair_nizkpop_Frodo640
#define crypto_nizkpop_verify         crypto_nizkpop_verify_Frodo640

#include "kem.c"
#include "noise.c"
#include "zkpop.c"
#if defined(USE_REFERENCE)
#include "frodo_macrify_reference.c"
#else
#include "frodo_macrify.c"
#endif
