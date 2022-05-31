/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: setting parameters to test FrodoKEM-1344
*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../src/api_frodo1344.h"


#define SYSTEM_NAME    "FrodoKEM-1344"

#define crypto_kem_keypair            crypto_kem_keypair_Frodo1344
#define crypto_kem_enc                crypto_kem_enc_Frodo1344
#define crypto_kem_dec                crypto_kem_dec_Frodo1344
#define crypto_kem_keypair_nizkpop    crypto_kem_keypair_nizkpop_Frodo1344
#define crypto_nizkpop_verify         crypto_nizkpop_verify_Frodo1344
#define shake                         shake256

#include "test_kem.c"
