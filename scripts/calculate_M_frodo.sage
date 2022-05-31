# Frodo parameter sets (n, q, beta, security strength in bits)
parameters = [
  (640, 2**15, 12, 128),
  (976, 2**16, 10, 192),
  (1344, 2**16, 6, 256)
]

for n,q,beta,secbits in parameters:
  # find biggest possible k
  for k in range(n):
    logprob = (log(1-2**-n) + log(2*beta+1)*(n-k) + log(q)*(k-n) + log(sum([ceil((2*beta+1)/(2**t))**n for t in range(int(log(q)/log(2)))])))/log(2)
    if -logprob < secbits:
      break
  print(f"Frodo-{n} ({secbits}-bit security)")
  print(f"k = {k}")
  sigma = n*16 # nbar is 8 for all parameter sets
  print(f"sigma = {sigma}")
  # find smallest possible M 
  M = sigma
  while True:
    prob = binomial(M-k, M-sigma) / binomial(M, M-sigma)
    if prob < 2**-secbits:
      break
    M += 1
  print(f"M = {M}")
  print("")
