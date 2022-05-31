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
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 1000

uint64_t t[NTESTS];
uint8_t seed[KYBER_SYMBYTES] = {0};

int main()
{
  unsigned int i;
  int rv;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t *zkpop;
  size_t zkpop_size;
  
  for(i=0;i<NTESTS;i++) {
    if ((rv = crypto_kem_keypair_zkpop(pk, sk, &zkpop, &zkpop_size)) == 0)
    {
      t[i] = zkpop_size;
      free(zkpop);
    } else {
      printf("bad kyber_keypair_zkpop: %d\n", rv);
      return -1;
    }
  }
  print_size_results("kyber_keypair_zkpop: ", t, NTESTS);
  
  return 0;
}
