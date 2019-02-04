#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G


#function decimalToVector
#input: numbers n and r (0 <= n<2**r)
#output: a string v of r bits representing n
def decimalToVector(n,r): 
    v = []
    for s in range(r):
        v.insert(0,n%2)
        n //= 2
    return v

#################################################################################

import numpy as np

# Functions for Hamming codes

def message(a):
    r = 2
    while 2**(r) -2*r -1 < len(a):
        r += 1
    k = 2**(r) -r -1
    print("k", k)
    print("r", r)
    out = []
    
    out += (decimalToVector(len(a), r))
    out += a
    while len(out) < k:
        out.append(0)
    return out

def hammingEncoder(m):
    r = 2
    while len(m) > 2**r-r-1:
        r += 1
    if len(m) < 2**r-r-1:
        return []
    m = np.array(m)
    hamming = hammingGeneratorMatrix(r)
    out = m.dot(hamming)
    out = out.tolist()
    return out

    





#################################################################################

# Functions for repetition codes

def repetitionEncoder(m,n):
    return m*n



def reperitionDecoder(v):
    total = sum(v)
    if total < len(v) // 2:
        return [0]
    elif total > len(v) // 2:
        return [1]
    else:
        return []


