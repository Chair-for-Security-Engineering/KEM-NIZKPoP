#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kem.h"
#include "randombytes.h"
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "zkpop.h"

#define NTESTS 100

uint64_t t[NTESTS];
uint8_t seed[KYBER_SYMBYTES] = {0};

int main()
{
  unsigned int i;
  int rv;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
  uint8_t key_a[CRYPTO_BYTES];
  uint8_t key_b[CRYPTO_BYTES];
  uint8_t *zkpop;
  size_t zkpop_size;
  
  printf("testing zkpop\n");
  for (i = 0; i < NTESTS; i++)
  {
    rv = crypto_kem_keypair_nizkpop(pk, sk, &zkpop, &zkpop_size);
    if (rv != 0)
    {
      printf("bad zkpop gen: %d\n", rv);
      return -1;
    }
    
    rv = crypto_nizkpop_verify(pk, zkpop, zkpop_size);
    if (rv != 0)
    {
      printf("bad zkpop verify: %d\n", rv);
      return -1;
    }
    
    uint8_t buf[4];
    uint32_t *pos = (uint32_t*)buf;
    randombytes(buf, 4);
    uint8_t mask = (1<<(*pos % 8));
    *pos >>= 3;
    zkpop[*pos % zkpop_size] ^= mask;
    
    rv = crypto_nizkpop_verify(pk, zkpop, zkpop_size);
    if (rv == 0)
    {
      printf("bad zkpop verify, expected for random bit flip at %u, %02x: %d\n", *pos % zkpop_size, mask, rv);
      return -1;
    }
    
    //Bob derives a secret key and creates a response
    crypto_kem_enc(ct, key_b, pk);

    //Alice uses Bobs response to get her shared key
    crypto_kem_dec(key_a, ct, sk);

    if(memcmp(key_a, key_b, CRYPTO_BYTES)) {
      printf("ERROR keys\n");
      return 1;
    }
    
    randombytes(buf, 4);
    mask = (1<<(*pos % 8));
    *pos >>= 3;
    ct[*pos % CRYPTO_BYTES] ^= mask;
    
    //Alice uses Bobs response to get her shared key
    crypto_kem_dec(key_a, ct, sk);

    if(memcmp(key_a, key_b, CRYPTO_BYTES) == 0) {
      printf("ERROR keys: random bit flip should result in decaps error\n");
      return 1;
    }
    free(zkpop);
  }
  printf("good.\n");

  return 0;
}
