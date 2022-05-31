/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: benchmarking/testing KEM scheme
*********************************************************************************************/

#include "../src/random/random.h"
#include "../src/sha3/fips202.h"
#include "speed_print.h"
#include "cpucycles.h"

#ifdef DO_VALGRIND_CHECK
#include <valgrind/memcheck.h>
#endif

#ifdef DO_VALGRIND_CHECK
#define KEM_TEST_ITERATIONS   1
#else
#define KEM_TEST_ITERATIONS   1
#endif
#if ZKPOP_N >= 8192
#define KEM_BENCH_ITERATIONS  3
#else
#define KEM_BENCH_ITERATIONS  100
#endif
#define KEM_BENCH_SECONDS     1


static int kem_test(const char *named_parameters, int iterations) 
{
  int rc;
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t ss_encap[CRYPTO_BYTES], ss_decap[CRYPTO_BYTES];
    uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
    unsigned char bytes[4];
    uint32_t* pos = (uint32_t*)bytes;
    uint8_t Fin[CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES];
    uint8_t *zkpop;
    size_t zkpop_size;

    #ifdef DO_VALGRIND_CHECK
        if (!RUNNING_ON_VALGRIND) {
            fprintf(stderr, "This test can only usefully be run inside valgrind.\n");
            fprintf(stderr, "valgrind frodo640/test_KEM (or frodo976 or frodo1344)\n");
            exit(1);
        }
    #endif

    printf("\n");
    printf("=============================================================================================================================\n");
    printf("Testing correctness of key encapsulation mechanism (KEM), system %s, tests for %d iterations\n", named_parameters, iterations);
    printf("=============================================================================================================================\n");

    for (int i = 0; i < iterations; i++) {
        if (i < iterations/2) // test NIZKPoP
        {
          if ((rc = crypto_kem_keypair_nizkpop(pk, sk, &zkpop, &zkpop_size)) != 0)
          {
              printf("\n ERROR -- Key/ZKPoP generation failed: %d!\n", rc);
            return false; 
          }
          if ((rc = crypto_nizkpop_verify(pk, zkpop, zkpop_size)) != 0)
          {
              printf("\n ERROR -- ZKPoP verification failed: %d!\n", rc);
            return false; 
          }
          // testing zkpop verification after changing random bits
          randombytes(bytes, 4);
          *pos %= 16000;
          zkpop[*pos] ^= 0xFF;
          if (crypto_nizkpop_verify(pk, zkpop, zkpop_size) == 0)
          {
              printf("\n ERROR -- ZKPoP verification should fail for random bit flip at position %d (size %lu)!\n", *pos, zkpop_size);
              printf("%02x%02x%02x\n", zkpop[*pos-1], zkpop[*pos], zkpop[*pos+1]);
            return false; 
          }
          free(zkpop);
        } else { // test usual keygen
          crypto_kem_keypair(pk, sk);
        }
        
        crypto_kem_enc(ct, ss_encap, pk);
        crypto_kem_dec(ss_decap, ct, sk);
#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_DEFINED(ss_encap, CRYPTO_BYTES);
        VALGRIND_MAKE_MEM_DEFINED(ss_decap, CRYPTO_BYTES);
#endif
        if (memcmp(ss_encap, ss_decap, CRYPTO_BYTES) != 0) {
            printf("\n ERROR -- encapsulation/decapsulation mechanism failed!\n");
	        return false; 
        }
        
        // Testing decapsulation after changing random bits of a random 16-bit digit of ct
        randombytes(bytes, 4);
        *pos %= CRYPTO_CIPHERTEXTBYTES/2;
        if (*pos == 0) {
            *pos = 1;
        }
        ((uint16_t*)ct)[*pos] ^= *pos;
        crypto_kem_dec(ss_decap, ct, sk);
#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_DEFINED(ss_decap, CRYPTO_BYTES);
#endif

        // Compute ss = F(ct || s) with modified ct
        memcpy(Fin, ct, CRYPTO_CIPHERTEXTBYTES);
        memcpy(&Fin[CRYPTO_CIPHERTEXTBYTES], sk, CRYPTO_BYTES);
        shake(ss_encap, CRYPTO_BYTES, Fin, CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);
        
#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_DEFINED(ss_encap, CRYPTO_BYTES);
#endif
        if (memcmp(ss_encap, ss_decap, CRYPTO_BYTES) != 0) {
            printf("\n ERROR -- changing random bits of the ciphertext should cause a failure!\n");
	        return false;
        }
    }
    printf("Tests PASSED. All session keys matched.\n");
    printf("\n\n");

    return true;
}

static void kem_bench(void) 
{
    size_t i;
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t *zkpop;
    size_t zkpop_size;
    int rv;

    uint64_t s[KEM_BENCH_ITERATIONS];
    uint64_t t[KEM_BENCH_ITERATIONS];
    
    for (i = 0; i < KEM_BENCH_ITERATIONS; i++)
    {
      t[i] = cpucycles();
      if ((rv = crypto_kem_keypair_nizkpop(pk, sk, &zkpop, &zkpop_size)) != 0)
      {
        printf("bad crypto_kem_keypair_nizkpop gen: %d\n", rv);
        return;
      }
      free(zkpop);
    }
    print_results("crypto_kem_keypair_nizkpop: ", t, KEM_BENCH_ITERATIONS);
    
#if ZKPOP_N < 8192
    for(i = 0; i < KEM_BENCH_ITERATIONS; i++) 
    {
      t[i] = cpucycles();
      if ((rv = crypto_kem_keypair_nizkpop(pk, sk, &zkpop, &zkpop_size)) != 0)
      {
        printf("bad crypto_kem_keypair_nizkpop gen: %d\n", rv);
        return;
      }
      s[i] = zkpop_size;
      if ((rv = crypto_nizkpop_verify(pk, zkpop, zkpop_size)) == 0)
      {
        free(zkpop);
      } else {
        printf("bad crypto_nizkpop_verify: %d\n", rv);
        return;
      }
    }
    print_results("zkpop gen plus verify: ", t, KEM_BENCH_ITERATIONS);
    
    print_size_results("proof: ", s, KEM_BENCH_ITERATIONS);
#endif
}


int main(int argc, char **argv) 
{
    int OK = true;

    OK = kem_test(SYSTEM_NAME, KEM_TEST_ITERATIONS);
    if (OK != true) {
        goto exit;
    }

    if ((argc > 1) && (strcmp("nobench", argv[1]) == 0)) {}
    else {
        kem_bench();
    }

exit:
    return (OK == true) ? EXIT_SUCCESS : EXIT_FAILURE;
}
