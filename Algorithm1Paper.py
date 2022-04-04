# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 18:21:27 2021

@author: alexa
"""
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy import optimize  as op


def f(x):
    '''
    
    Objective function
    f(x)=(1/2)*||Ax-b||^2

    Parameters
    ----------
    x : vector 
        Argument .

    Returns
    -------
    value of x.

    '''
    
    return 0.5*np.linalg.norm(np.matmul(A,x)-b)**2
    
def gradf(x):
    '''
    

    Parameters
    ----------
    x : numpy array
        Vector .

    Returns
    -------
    Jacobian of f at x Df(x), since f is a scalar field the Jacobian is a vector.

    '''   
    return (np.matmul(A,x)-b).dot(A)


def outer_Frank_Wolfe(steps,V0):
    '''
    

    Parameters
    ----------
    steps : int
        Number of steps we will do .
    stept :float 
        current step of the main Algorithm .
    V : numpy array
        List containing vertices of the decomposition
        into k vertices, row wise.
    x : numpy array 
        Current Iterate of the Main Alg.
    Beta : float
        smoothnes Parameter.
    lamba : numpy array 
        Contains the scalars of the decomposition into vertices.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''

    def Oracle(x,Aub,bub):
        '''
        A Linear Oracle for the subproblem 
    
        Parameters
        ----------
        tauk : numpy array 
            current iterate of subproblem.
        V : numpy array 
            vertices of decomposition.
        x : numpy array
            current iterate in main Alg.
        Beta : Float
            Smoothnes Parameter.
        stept : Float
            Current stepsize in main Alg .
        Aub : numpy array
            constraints for domain of Tau.
        Aeq : numpy array
            contraints for sum over components.
        bub : numpy array
            correspondent boundry to Aub.
        beq : TYPE
            correspondent boundry to Aeq.
    
        Returns
        -------
        numpy array
            new direction.
    
        '''
          
        c=gradf(x)
      
        return op.linprog(c,A_ub=Aub,b_ub=bub).x
    
    step2=1
    xk=V0

    Iterations=[xk]
    for i in range(steps):
        print(i,'FW')
        Aub=np.identity(d)
        # Aeq=np.zeros((k,k))
        # Aeq[0]=np.ones(k)
        bub=np.ones(d)
        #beq=np.zeros(k)
       # beq[0]=1
        w=Oracle(xk,Aub,bub)
        
        xk=(1-step2)*xk+step2*w

        Iterations.append(xk)
        step2=2/(i+2)
    return Iterations

def Algorithm1(x0,sizes,Beta):
    
    def NEP_Oracle1(x0,Beta,Sizesi):
        
        lamba=Beta*Sizesi
     
        G=x0-1/(2*lamba)*gradf(x0)
       
        v=np.zeros(len(x0))
        for i in range(len(G)):
            if G[i]<1/2:
                v[i]=0
            else:
                v[i]=1
        return v
                
        
        
        
        
    def adaptive(xs,steps,k,w):
        valnew=f(xs[-1] + steps[i]*(w-xs[-1]))
        valold=f(xs[-1])
        if valold<valnew:
            return 0
        else: 
            return steps[k]
    Iterations=[x0]
    for i in range(len(sizes)):
        print(i,'NEP')
        w=NEP_Oracle1(x0,Beta,sizes[i])
        
        step=adaptive(Iterations,sizes,i,w)
        x0=(1-step)*x0+step*w
        Iterations.append(x0)
    return Iterations


##############################################################################
##############################################################################
################### WE NOW SET THE STAGE #####################################
##############################################################################
##############################################################################
m,d=100,200
# Vertices =np.array( list(itertools.product([1,0], repeat=d)))
A=np.random.normal( size=(m, d))
#A=np.identity(d)
X=np.zeros(d)
Randvert=np.ones(d)
X=Randvert
X[:5]=0.5*np.ones(5)
b=np.matmul(A,X)
V0=np.ones(d)
V0[:6]=np.zeros(6)
print('Step1')
How_many=100
Sizes=[1/(i+2) for i in range(How_many)]
Numbers=[i for i in range(How_many)]

Beta=np.linalg.eigh(np.matmul(A,A.T))[0][-1]
###############################################################3
# Result1=Frank_Wolfe_NEP(V0,Vertices,Beta, Sizes)
# Iterations=Result1[1]
# Values=[f(Iterations[i]) for i in range(How_many)]
# plt.plot(Numbers,Values,label='Algorithm 2')
####################################################################
plt.ylim([0, 100])

Result2=np.array(outer_Frank_Wolfe(How_many,V0))
plt.xlabel('Iteration')
plt.ylabel('function value')
Values2=[f(Result2[i]) for i in range(How_many)]
plt.plot(Numbers,Values2,label='FWolfe')
########################################################################
Result3=np.array(Algorithm1(V0, Sizes,Beta))
Values3=[f(Result3[i]) for i in range(How_many)]
plt.plot(Numbers,Values3,label='Algorithm1')

plt.legend()
