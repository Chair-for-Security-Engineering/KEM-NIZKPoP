from math import log,ceil,floor
from natsort import natsorted
from tabulate import tabulate
import glob

# Frodo parameter sets (n, q, beta, security strength in bits, sigma, M)
parameters = [
  (640, 2**15, 12, 128, 10240, 13233, 9616),
  (976, 2**16, 10, 192, 15616, 19485, 15632),
  (1344, 2**16, 6, 256, 21504, 25986, 21520)
]

def idealsize(secpar, tau, N, M, sigma, beta, q):
  return 6*secpar + tau * (2*secpar + secpar*ceil(log(N,2)) + ceil(log(q,2)) * M) + ceil(log(2*beta+1,2)) * (M - sigma)

for n, q, beta, secbits, sigma, M, pksize in parameters:
  fnl = natsorted(glob.glob(f"frodo{n}_*"))
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
      res.append((I, median_keypair_nizkpop, median_total-median_keypair_nizkpop, (min_size)/1000, (median_size)/1000, idealsize(secbits, tau, N, M, sigma, beta, q)/8000,  (max_size)/1000, f"{{{N,tau}}}"))
    except TypeError:
      pass
  
  with open(f"frodocycles{n}", "w") as f:
    f.write(tabulate(res, headers=["index", "gencycles", "vercycles", "minsizemeasured", "mediansize", "idealsize", "maxsizemeasured", "label"], tablefmt="plain"))
    f.write("\n")
