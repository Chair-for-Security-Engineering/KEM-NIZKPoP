#ifndef ZKPOP__H
#define ZKPOP__H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#define LAZY 9

typedef ALIGNED_INT16(ZKPOP_M) bundle_t;
typedef ALIGNED_INT16(ZKPOP_M-ZKPOP_SIGMA) open_bundle_t;

#if ZKPOP_N > 65536
#error "Too many parties."
#elif ZKPOP_N > 256
typedef uint16_t party_t;
#else 
typedef uint8_t party_t;
#endif

#define crypto_kem_keypair_nizkpop KYBER_NAMESPACE(crypto_kem_keypair_nizkpop)
int crypto_kem_keypair_nizkpop(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES],
                    uint8_t **zkpop,
                    size_t *zkpop_size);

#define crypto_nizkpop_verify KYBER_NAMESPACE(crypto_nizkpop_verify)
int crypto_nizkpop_verify(const unsigned char *pk, const unsigned char *zkpop, unsigned long zkpop_size);


void generate_S_E(polyvec *S, polyvec *E, const bundle_t *vk, const uint8_t audit_bundles[(ZKPOP_M+7)/8]);
void generate_S_E_openbundle(polyvec *S, polyvec *E, open_bundle_t *openbundle, const bundle_t *vk, const uint8_t audit_bundles[(ZKPOP_M+7)/8]);
void generate_seeds(uint8_t *seedtree);
void sample_audit_parties (party_t audit_parties[ZKPOP_TAU], const uint8_t hB[KYBER_SYMBYTES]);
void sample_audit_bundles (uint8_t audit_bundles[(ZKPOP_M+7)/8], const uint8_t h[KYBER_SYMBYTES]);

#endif
