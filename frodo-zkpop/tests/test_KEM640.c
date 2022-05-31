/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: setting parameters to test FrodoKEM-640
*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../src/api_frodo640.h"


#define SYSTEM_NAME    "FrodoKEM-640"

#define crypto_kem_keypair            crypto_kem_keypair_Frodo640
#define crypto_kem_enc                crypto_kem_enc_Frodo640
#define crypto_kem_dec                crypto_kem_dec_Frodo640
#define crypto_kem_keypair_nizkpop    crypto_kem_keypair_nizkpop_Frodo640
#define crypto_nizkpop_verify         crypto_nizkpop_verify_Frodo640
#define shake                         shake128

#include "test_kem.c"
