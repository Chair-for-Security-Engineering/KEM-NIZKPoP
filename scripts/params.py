from scipy.stats import binom
import math

# Script to find parameters and corresponding proof sizes
# The sizes here *exclude* the public key

def sigSize(kappa, q_bits, beta_bits, sigma, M, N, tau):
    # Signature size formula from section 4.2
    sigSize = 0.0
    sigSize += 2*kappa  # salt
    sigSize += 4*kappa  # h_1, h_2

    proofSize = 0.0 # this will be the size of data sent per-repetition
    proofSize += 2*kappa  # commitment to seed of unopened party
    proofSize += kappa * math.ceil(math.log(N, 2))  # reveallist of seed tree
    proofSize += q_bits * M 

    sigSize += beta_bits * (M - sigma)
    sigSize += tau * proofSize
    return to_bytes(sigSize)


def to_bytes(x):
    return math.ceil(x/8.0)

def soundness(kappa, q_bits, sigma, M, N, tau):
    p1 = 2**kappa
    p2 = N**tau
    return(p2)

def print_helper(kappa, q_bits, beta_bits, sigma, M, N, tau, cost):
    stars = ''
    if abs(N-256) < 10:
        stars = '***'
    print("N={}, tau={}, q_bits={}, beta_bits={}, M={}, sigma={},  2^{:.3f} attack cost, SigSize: {:,}   {}".format(
        N, tau, q_bits, beta_bits, M, sigma, float(math.log(cost, 2)), sigSize(kappa, q_bits, beta_bits, sigma, M, N, tau), stars))


def paramsByParties(kappa, q_bits, beta_bits, M, sigma, parties):
    # Finds the parameters for a set of chosen N
    print("{}-bit security, q=2^{}, beta=2^{}, M = {}, sigma = {}".format(kappa, q_bits, beta_bits, M, sigma))
    for N in parties:
        for tau in range(2*kappa):
            cost, tau_1 = soundness(kappa, q_bits, sigma, M, N, tau)
            if cost >= 2**kappa:
                print_helper(kappa, q_bits, beta_bits, sigma, M, N, tau, cost)
                break


def paramsByRepetitions(kappa, q_bits, beta_bits, M, sigma, rounds):
    # Finds the parameters for a fixed set of chosen \tau
    print("{}-bit security, q=2^{}, beta=2^{}, M = {}, sigma = {}".format(kappa, q_bits, beta_bits, M, sigma))
    MAX_N = 2**16
    for tau in rounds:
        cost = 0
        for N in range(2, MAX_N):   # Find minimal N for this tau, giving 2**kappa security
            cost = soundness(kappa, q_bits, sigma, M, N, tau)
            if cost >= 2 ** kappa:
                break
            N = N+1
        if N == MAX_N and cost < 2**kappa:
            print("Unrealistic number of parties needed, tau = ", tau, " cost = ", math.log(cost,2))

        print_helper(kappa, q_bits, beta_bits, sigma, M, N, tau, cost)

print("-"*80)
logParties = list(range(6, 10))
parties = [int(2**n) for n in logParties]


def L1frodo640():
    print("L1: Frodo 640")
    sigma = 10240
    M =  sigma + 2993
    q_bits = 15
    beta_bits = 5

    print("-"*80)
    paramsByRepetitions(128, q_bits, beta_bits, M, sigma, range(8,90))
    print("-"*80)

def L1kyber512():
    print("L1: Kyber-512")
    sigma = 1024 
    M =  1280 
    q = 3329
    q_bits = math.ceil(math.log(q,2))
    beta_bits = 2

    print("-"*80)
    paramsByRepetitions(128, q_bits, beta_bits, M, sigma, range(8,50))
    print("-"*80)


def Regime1():
    print("Parameter Regime 1")
    sigma = 2048
    M =  2386 
    q = 2**32
    q_bits = math.ceil(math.log(q,2))
    beta_bits = math.log(3,2)

#    paramsByParties(128, q_bits, beta_bits, M, sigma, parties)
    print("-"*80)
    paramsByRepetitions(128, q_bits, beta_bits, M, sigma, range(8,50))
    print("-"*80)

def Regime2():
    print("Parameter Regime 2")
    sigma = 4096
    M =  4419
    q = 2**61
    q_bits = math.ceil(math.log(q,2))
    beta_bits = math.log(3,2)

#    paramsByParties(128, q_bits, beta_bits, M, sigma, parties)
    print("-"*80)
    paramsByRepetitions(128, q_bits, beta_bits, M, sigma, range(8,50))
    print("-"*80)

L1frodo640()
L1kyber512()
Regime1()
Regime2()
