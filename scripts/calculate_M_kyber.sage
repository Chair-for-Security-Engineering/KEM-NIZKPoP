# Kyber parameter sets (n, l, q, eta1, security strength in bits)
parameters = [
  (256, 2, 3329, 3, 128),
  (256, 3, 3329, 2, 192),
  (256, 4, 3329, 2, 256)
]

for n, l, q, eta1, secbits in parameters:
  # find biggest possible k
  for k in range(n*l):
    logprob = log(2*eta1+1,2) * l*n + log((2*eta1+1)/q,2)*(l*n-k)
    if -logprob < secbits:
      break
  print(f"Kyber-{n*l} ({secbits}-bit security)")
  print(f"k = {k}")
  sigma = n*l*2
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
