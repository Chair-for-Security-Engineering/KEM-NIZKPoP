from math import log,ceil,floor
from natsort import natsorted
from tabulate import tabulate
import glob

# Kyber parameter sets (n, l, q, eta1, security strength in bits, sigma, M, pksize)
parameters = [
  (256, 2, 3329, 3, 128, 1024, 1280),
  (256, 3, 3329, 2, 192, 1536, 1870),
  (256, 4, 3329, 2, 256, 2048, 2493)
]

def idealsize(secpar, tau, N, M, sigma, eta, q):
  return 6*secpar + tau * (2*secpar + secpar*ceil(log(N,2)) + ceil(log(q,2)) * M) + ceil(log(2*eta+1,2)) * (M - sigma)
  

for n, l, q, eta1, secbits, sigma, M in parameters:
  pksize = 384*l + 32
  fnl = natsorted(glob.glob(f"kyber{n*l}*"))
  res = []
  keygens = []
  for I,fn in enumerate(fnl):
    tmp = fn.split("_")
    N = int(tmp[1][1:])
    tau = int(tmp[2][3:])
    with open(fn, "r") as f:
      c = f.readlines()
    
    median_keypair_nizkpop = None
    for i,line in enumerate(c):
      if line.find("crypto_kem_keypair_nizkpop") != -1:
        median_keypair_nizkpop = int(c[i+1].split()[1])
    
    median_total = None
    for i,line in enumerate(c):
      if line.find("zkpop gen plus verify") != -1:
        median_total = int(c[i+1].split()[1])
        
    median_size = None
    min_size = None
    max_size = None
    for i,line in enumerate(c):
      if line.find("proof:") == 0:
        min_size = int(c[i+1].split()[1])
        median_size = int(c[i+2].split()[1])
        max_size = int(c[i+4].split()[1])
    try:
      res.append((I, median_keypair_nizkpop, median_total-median_keypair_nizkpop, min_size/1000, (median_size)/1000, idealsize(secbits, tau, N, M, sigma, eta1, q)/8000, (max_size)/1000, f"{{{N,tau}}}"))
    except TypeError:
      pass
  
  with open(f"kybercycles{n*l}", "w") as f:
    f.write(tabulate(res, headers=["index", "gencycles", "vercycles", "minsizemeasured", "mediansize", "idealsize", "maxsizemeasured", "label"], tablefmt="plain"))
    f.write("\n")
