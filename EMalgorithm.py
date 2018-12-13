# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 17:51:01 2017

@author: victor
"""

import time, math, random
import numpy as np

N=190   #number of stations
D=30	#number of days
T=48    #number of time window
clusters=3    #number of desired clusters
L=2     #day cluster
filelocation= "EMdata30_1.txt"
MaxIter=100 # maximum runs of iterations


def poisson(k,mean):
        return (np.exp(-mean))*(mean)**k/math.factorial(k)

BICplot=[]
AICplot=[]
loglikelihoodplot=[]

for K in range(1,clusters):
    print 'cluster:',K
    W=np.zeros((D,2),dtype=np.int)
    for d in range(0,D,7):    
        W[d:d+5,0]=1
        W[(d+5):(d+7),1]=1  
    
    inFile = open(filelocation)
    tmpIn = inFile.readline().strip().split("\t")
    X=np.zeros((N,D,T),dtype=np.int)
    n=0
    d=0
    for x in inFile:
        tmpIn = x.strip().split("\t")
        X[n][d][0:T]=tmpIn[2:2+T]
        d=d+1
        if d==D: 
            d=0
            n=n+1
        else: continue
    
    #initialization of alpha
    alpha=[]
    for s1 in range(N):
        SUM=np.sum(X[s1,n1,t1] for n1 in range(D) for t1 in range(T))
        alpha.append(float(SUM)/D/T)
    
    #initialization of pi
    pi={}
    currentpi=[]
    for p in range(K):
        #currentpi=[0.33,0.21,0.19,0.17,0.10]
        currentpi.append(1.0/K)
    pi[0]=currentpi
    
    #initialization of lambda
    Dweekdays=np.sum(W[:,0])
    Dweekends=np.sum(W[:,1])
    lam={}
    lam[0]=np.full((K,L,T),float(D*T)/(T*(Dweekdays+Dweekends)))
    for k in range (K):
        for l in range (L):
            for t in range (T/2):
                r=random.random()
                lam[0][k,l,t]=lam[0][k,l,t]+r
                lam[0][k,l,t+T/2]=lam[0][k,l,t+T/2]-r

    t_sk=np.zeros((N,K))
    diff=10
    precision=0.00001
    loglikelihoodlist=[]
    iter=0
    
    while diff>=precision and iter<=MaxIter:
    #E step 
        for s2 in range(N):
            lognumerator=[]
            for k2 in range(K):
                term=[]
                for l2 in range(L):
                    for d2 in range(D):
                        for t2 in range (T):
                                term.append(np.log(poisson(X[s2,d2,t2],alpha[s2]*lam[iter][k2,l2,t2])**W[d2,l2]))
                lognumerator.append(np.sum(term)+np.log(pi[iter][k2]))
            MAX=max(lognumerator)
            logdenominator=[x-MAX for x in lognumerator]
            for k3 in range(K):
                t_sk[s2,k3]=np.exp(logdenominator[k3])/np.sum(np.exp(logdenominator[i]) for i in range(K))
    #M step
        #update pi
        currentpi=[]
        for k4 in range(K):
            currentpi.append(np.nansum(t_sk[:,k4])/N)
        pi[iter+1]=currentpi
    
        #update lambda
        lam[iter+1]=np.ones((K,L,T))# 
        for k5 in range(K):
            for l5 in range (L):
                for t5 in range(T):
                    lam[iter+1][k5,l5,t5]=np.nan_to_num(np.sum(t_sk[s5,k5]*float(W[d5,l5])*X[s5,d5,t5] for s5 in range(N) for d5 in range(D))/(np.sum((t_sk[s5,k5])*alpha[s5] for s5 in range (N))*float(np.sum(W[:,l5]))))
    
    #calculate loglikelihood
        loglikelihoodsum=[]
        for s2 in range(N):
            lognumerator=[]
            for k2 in range(K):
                term=[]
                for l2 in range(L):
                    for d2 in range(D):
                        for t2 in range (T):
                                term.append(np.log(poisson(X[s2,d2,t2],alpha[s2]*lam[iter+1][k2,l2,t2])**W[d2,l2]))
                lognumerator.append(np.sum(term)+np.log(pi[iter+1][k2]))
            loglikelihoodsum.append(np.log(np.sum(np.exp(lognumerator[i]) for i in range(K))))        
        nonzerologlikelihoodsum=[e  for i,e in enumerate(loglikelihoodsum) if e!=-np.inf]
        loglikelihood=np.sum(nonzerologlikelihoodsum[ss] for ss in range(len(nonzerologlikelihoodsum)))
        loglikelihoodlist.append([iter,loglikelihood])
        print 'finish calculating loglikelihood'
        
        if iter==0:
            diff=10
        else:
            diff=np.abs(loglikelihoodlist[iter][1]-loglikelihoodlist[iter-1][1])
        iter=iter+1
        
    #plot convergence figure 
    BICvalue=-2*loglikelihoodlist[-1][1]+3*K*np.log(T)
    BICplot.append([K,BICvalue])
    AICvalue=-2*loglikelihoodlist[-1][1]+2*3*K
    AICplot.append([K,AICvalue])
    loglikelihoodplot.append([K,loglikelihoodlist[-1][1]])
 







