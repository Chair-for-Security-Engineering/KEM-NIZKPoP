#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "kem.h"
#include "kex.h"
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "cpucycles.h"
#include "speed_print.h"
#include "zkpop.h"
#include "cbd.h"

#if ZKPOP_N >= 8192
#define NTESTS 3
#else
#define NTESTS 100
#endif

uint64_t t[NTESTS];
uint64_t s[NTESTS];
uint8_t seed[KYBER_SYMBYTES] = {0};

/* Dummy randombytes for speed tests that simulates a fast randombytes implementation
 * as in SUPERCOP so that we get comparable cycle counts */
void randombytes(__attribute__((unused)) uint8_t *r, __attribute__((unused)) size_t len) {
  return;
}

int main()
{
  unsigned int i;
  int rv;
  uint8_t *zkpop;
  size_t zkpop_size;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  
#ifndef PROFILE
  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    crypto_kem_keypair(pk, sk);
  }
  print_results("kyber_keypair: ", t, NTESTS);
  
  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    if ((rv = crypto_kem_keypair_nizkpop(pk, sk, &zkpop, &zkpop_size)) == 0)
    {
      free(zkpop);
    } else {
      printf("bad crypto_kem_keypair_nizkpop: %d\n", rv);
      return -1;
    }
  }
  print_results("crypto_kem_keypair_nizkpop: ", t, NTESTS);
#endif

  for(i=0;i<NTESTS;i++) {
    t[i] = cpucycles();
    if ((rv = crypto_kem_keypair_nizkpop(pk, sk, &zkpop, &zkpop_size)) != 0)
    {
      printf("bad crypto_kem_keypair_nizkpop gen: %d\n", rv);
      return -1;
    }
    s[i] = zkpop_size;
    if ((rv = crypto_nizkpop_verify(pk, zkpop, zkpop_size)) == 0)
    {
      free(zkpop);
    } else {
      printf("bad crypto_nizkpop_verify: %d\n", rv);
      return -1;
    }
  }
  print_results("zkpop gen plus verify: ", t, NTESTS);
  
  print_size_results("proof: ", s, NTESTS);

  return 0;
}
