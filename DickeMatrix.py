#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 13:33:57 2021

@author: michel
"""
import numpy as np
import matplotlib.pyplot as p
import time

def basisDimension(N, gamma):
    dim = int((gamma + 1)*(N + 2)*(N + 1)/2)    # dimension del espacio de Hilbert
    return dim

def basisVectors(N, gamma):
    nBasis = []                                 # base de partículas
    
    for i in range(3):                          # N partículas en un único nivel
        vec = np.zeros(3)
        vec[i] = N
        nBasis.append(vec.copy())
        
    for i in range(1,N):                        # N = N1 + N2
        N1 = N - i
        N2 = i
        for j in range(3):                      # un nivel vacío
            vec = np.zeros(3)
            vec[j] = N1
            vec[j-1] = N2
            nBasis.append(vec.copy())
        if N1 != 1:                             # si N1 puede descomponerse en suma
            for j in range(1,N1):               # N1 = N11 + N12
                N11 = N1 - j
                N12 = j
                vec = np.zeros(3)
                vec[0] = N11
                vec[1] = N12
                vec[2] = N2
                nBasis.append(vec.copy())
    basis = []
    
    for i in range(len(nBasis)):
        for j in range(gamma+1):
            gvec = np.append(nBasis[i], j)
            basis.append(gvec.copy())
    
    return basis

def Hact(vi, vj, wf, wn, D, f, N, E):
    H = 0
    if vi[0] == vj[0] and vi[1] == vj[1] and vi[2] == vj[2] and vi[3] == vj[3]:
        H += wf*vi[3]
        
        for alpha in range(3):
            H += wn[alpha]*vi[alpha]

    if vi[0] == vj[0] and vi[1] == vj[1] and vi[2] == vj[2]:
        if (vi[3] - 2) == vj[3]:
            H += D*np.sqrt(vi[3]*(vi[3] - 1))
        elif (vi[3] + 2) == vj[3]:
            H += D*np.sqrt((vi[3] + 1)*(vi[3] + 2))
        elif vi[3] == vj[3]:
            H += D*(2*vi[3] + 1)
            
        elif (vi[3] - 1) == vj[3]:
            H += E*np.sqrt(vi[3]/N)
        elif (vi[3] + 1) == vj[3]:
            H += E*np.sqrt((vi[3] + 1)/N)

    c = [2,1,0]
    Om = np.zeros((3,3))
    for b in [1,2]:
        for a in range(b):
            Om[a][b] = np.sqrt((wn[b] - wn[a])*f[a+b-1]*D)
            if vi[c[a+b-1]] == vj[c[a+b-1]]:
                if (vi[3] - 1) == vj[3]:
                    if (vi[a] + 1) == vj[a] and (vi[b] - 1) == vj[b]:
                        H += Om[a][b]*np.sqrt((vi[a] + 1)*vi[b]*vi[3]/N)
                    elif (vi[a] - 1) == vj[a] and (vi[b] + 1) == vj[b]:
                        H += Om[a][b]*np.sqrt((vi[b] + 1)*vi[a]*vi[3]/N)
                elif (vi[3] + 1) == vj[3]:
                    if (vi[a] + 1) == vj[a] and (vi[b] - 1) == vj[b]:
                        H += Om[a][b]*np.sqrt((vi[a] + 1)*vi[b]*(vi[3] + 1)/N)
                    elif (vi[a] - 1) == vj[a] and (vi[b] + 1) == vj[b]:
                        H += Om[a][b]*np.sqrt((vi[b] + 1)*vi[a]*(vi[3] + 1)/N)
    
    return H

def H(gamma, basis, wf, wn, D, f, N, E):
    dim = basisDimension(N, gamma)
    Hmatrix = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            Hmatrix[i][j] = Hact(basis[i], basis[j], wf, wn, D, f, N, E)
    return Hmatrix

def eigenvalues(Hmatrix):
    eigvals = np.linalg.eigvals(Hmatrix)
    return eigvals

def eigvalHistogram(eigvals, n):
    
    p.hist(eigvals, density = True, bins = n, label = 'eigenvals histogram')
    p.show()

def rutine(Hmatrix):
    sol = []
    for i in range(len(Hmatrix)):
        for j in range(len(Hmatrix)):
            sol.append(Hmatrix[i][j] == Hmatrix[j][i])
    return sol

t0 = time.time()

N = 12
gamma = 1
wf = 1
wn = [1, 1.1, 2.1]
D = 3
f = [6, 6, 6]
E = 0
basis = basisVectors(N, gamma)
Hmatrix = H(gamma, basis, wf, wn, D, f, N, E)
#print(basis)
print(len(basis))
print(len(Hmatrix))
print(Hmatrix[181][181])
eigvals = eigenvalues(Hmatrix)
#print(eigvals)
eigvalHistogram(eigvals, 100)
#sol = rutine(Hmatrix)
#print(sol)
d = time.time() - t0
print("Timelapse: ", d, " seconds")


'''def prodEsc(vi, vj):
    if vi[0] == vj[0] and vi[1] == vj[1] and vi[2] == vj[2] and vi[3] == vj[3]:
        prod = 1
    else:
        prod = 0
    return prod

def OP_a(vec):
    if vec[3] != 0:
        veca = vec.copy()
        prod = np.sqrt(veca[3])
        veca[3] -= 1
    else:
        prod = 0
        veca = np.zeros(4)
    return prod, veca

def OP_ad(vec):
    veca = vec.copy()
    prod = np.sqrt(veca[3] + 1)
    veca[3] += 1
    return prod, veca

def OP_Sigma(vec, alpha, beta):
    veca = vec.copy()
    if alpha == beta:
        prod = veca[alpha]
    else:
        if veca[beta] != 0:
            prod = np.sqrt((veca[alpha] + 1)*veca[beta])
            veca[alpha] += 1
            veca[beta] -= 1
        else:
            prod = 0
            veca = np.zeros(4)
    return prod, veca
    
def Hact(vi, vj, w_f, w_n, D, W_n, N, E):
    a, via = OP_a(vi)
    ad, viada = OP_ad(via)
    ifis = prodEsc(viada, vj)
    Wada = ifis*ad*a*w_f
    
    ###
    
    Wn = 0
    for alpha in range(3):
        nalpha, vialpha = OP_Sigma(vi, alpha, alpha)
        ifis = prodEsc(vialpha, vj)
        Wn += ifis*nalpha*w_n[alpha]

    ###
    
    a, via = OP_a(vi)
    a2, viaa = OP_a(via)
    ifis = prodEsc(viaa, vj)
    aa = ifis*a2*a
    
    ad, viad = OP_ad(vi)
    ad2, viadd = OP_ad(viad)
    ifis = prodEsc(viadd, vj)
    adad = ifis*ad2*ad
    
    a, via = OP_a(vi)
    ad, viada = OP_ad(via)
    ifis = prodEsc(viada, vj)
    ada = ifis*ad*a
    
    ad, viad = OP_ad(vi)
    a, viaad = OP_a(viad)
    ifis = prodEsc(viaad, vj)
    aad = ifis*a*ad
    
    Dada = (aa + adad + ada + aad)*D
    
    ###
    
    a, via = OP_a(vi)
    ad, viad = OP_ad(vi)
    Sn = 0
    for alpha in range(3):
        for beta in range(3):
            if alpha != beta:
                nab, viab = OP_Sigma(via, alpha, beta)
                ifis = prodEsc(viab, vj)
                Sn += ifis*nab
                
                nba, viba = OP_Sigma(via, beta, alpha)
                ifis = prodEsc(viba, vj)
                Sn += ifis*nba
                
                ndab, vidab = OP_Sigma(viad, alpha, beta)
                ifis = prodEsc(vidab, vj)
                Sn += ifis*ndab
                
                ndba, vidba = OP_Sigma(viad, beta, alpha)
                ifis = prodEsc(vidba, vj)
                Sn += ifis*ndba
                
                Sn = Sn*W_n[alpha, beta]/np.sqrt(N)
    
    ###
    
    a, via = OP_a(vi)
    ifis = prodEsc(via, vj)
    Ea = ifis*a
    ad, viad = OP_ad(vi)
    ifis = prodEsc(viad, vj)
    Ea += ifis*ad
    Ea = Ea*E/np.sqrt(N)
    
    Hij = Wada + Wn + Dada + Sn + Ea
    
    return Hij

def H(N, gamma, basis, w_f, w_n, D, W_n, E):
    dim = basisDimension(N, gamma)
    Hmatrix = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            Hmatrix[i][j] = Hact(basis[i], basis[j], w_f, w, D, W_n, N, E)
    return Hmatrix'''