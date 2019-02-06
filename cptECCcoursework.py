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

def vectorToDecimal(v):
    x = 1
    total = 0
    for i in range(len(v)-1, -1, -1):
        total += v[i]*x
        x = x * 2
    return total


#################################################################################

import numpy as np
import math

# Functions for Hamming codes

def message(a):
    r = 2
    while 2**(r) -2*r -1 < len(a):
        r += 1
    k = 2**(r) -r -1
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
    out = np.mod(out, np.full_like(out, 2))

    out = out.tolist()
    return out

def hammingDecoder(v):
    r = 2
    while len(v) > 2**r-1:
        r += 1
    if len(v) < 2**r-1:
        return []
    
    Htranspose = [decimalToVector(i, r) for i in range(1, 2**r)]
    v=np.array(v)
    Htranspose = np.array(Htranspose)

    
    errorpos = v.dot(Htranspose)
    errorpos = np.mod(errorpos, np.full_like(errorpos, 2))

    if not errorpos.any():
        return v.tolist()
    errorpos_decimal = vectorToDecimal(errorpos.tolist())-1

    v = v.tolist()
    v[errorpos_decimal] = (v[errorpos_decimal]+1)%2

    return v

    

    


def messageFromCodeword(c):
    r = 2
    while len(c) > 2**r-1:
        r += 1
    if len(c) < 2**r-1:
        return []
    output = []
    for i,n in enumerate(c):
        if not math.isclose(math.log2(i+1)%1, 0):
            output.append(n)
    return output


def dataFromMessage(m):
    r = 2
    while len(m) > 2**r-r-1:
        r += 1
    if len(m) < 2**r-r-1:
        return []
    length = vectorToDecimal(m[:r])

    if len(m) - r < length:
        return []
    message = m[r:r+length]
    return message






#################################################################################

# Functions for repetition codes

def repetitionEncoder(m,n):
    return m*n



def repetitionDecoder(v):
    total = sum(v)
    if total < len(v) / 2:
        return [0]
    elif total > len(v) / 2:
        return [1]
    else:
        return []


