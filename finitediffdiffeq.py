# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 16:03:05 2021
@author: MarlaJahari
"""
import numpy as np
import matplotlib.pyplot as ourplot

"PROBLEM: To solve the given diff eqn using finite difference method"
"y2=y'+x(x-4)  p(x)=1 q(x)=0 r(x)=x(x-4) 0<x<4 y(0)=y(4)=0"
'[a,b]=[0,4]'

a=0                   #upper limit
b=4                   #lower limit
n=int(input('no of steps YOU WISH'))
h=(b-a)/n
p=1                   #co-efficient of first derivative
q=0                   #co-efficient of y

'(2-p*h)y[i+1]+ (-4)y[i] + (2+p*h)y[i-1]=2rh^2 (substituting our equation in FDM)'

c=[]                 #vector to store the constants on right 

for i in range(1,n,1):
    r=2*i*h*((i*h)-4)#since x=0+i*h
    C= r*h**2        #the constant part of equation
    c.append(C)      #producing elements in loop and adding them to c

A=-4                 #(-2*q*h**2-4) co-efficient of y[i] throughout the loop   
B=2-p*h              #co-efficient of y[i+1] "            
E=2+p*h              #co-efficient of y[i-1] "
x=[]                 #for storing A for diff values of i
y=[]                 #for storing B for "
z=[]                 #for storing E for "             

for i in range(n-1): 
    x.append(A)      #adding B to vector x for diff i's

for i in range(n-2):
    y.append(B)      #adding A to vector y for diff i's
    z.append(E)

def tridiag(x,y,z, k1=-1,k2=0,k3=1): #creating a function for triD matrix
    
     return np.diag(x,k2) + np.diag(z,k1)+ np.diag(y,k3) #main+off+off

M=tridiag(x, y, z)                   #built matrix
P=np.array([[c]]).T                  #converting const array to soln column vector                      
print("Our co-efficient matrix:")
print(M)                             #printing matrix (*not printed with the command for spacing)
print("Our constant matrix:")
print(P)                             #printing soln vector


'Algorithm for Thomas method to solve tri-diagonal matrixes'

"""
here consider our reference matrix is of form
   
   [b1 c1 0
   a2 b2 c2
   0  a3 b2]

"""

v=len(c)                   #extracting the length of const mtrx
f=np.zeros(v-1)            #creating an array of zeros length v-1
g=np.zeros(v)                 
d=np.zeros(v)             

'~Matrix Operations~'

# for an eq b1x1+ c1x2=r1, div throughout by b1(1st element of main D)>>

f[0]=y[0]/x[0] #dividing 1st elm(Up-D) by 1st elm(M-D)
g[0]=c[0]/x[0]   
               

#row2- eqn a2x1+b2x2+c2x3, to eliminate a2x1>> 
#take the 1st element of off-D below(a2) and multiply by the 1st row
#repeat it down the rows until  all the MD elms=1, lD elms=0, uD elms=k
#we divide throughout by b2-a2*(k) k=previous uD modified elm>>

for i in range(1,v-1): 
              
      f[i]=y[i]/(x[i]-z[i-1]*f[i-1]) #co-efficients of xi throughout 

for i in range(1,v):
    
      #now we modify const Vectr in sync with the rest of our operations
      #modifying constant Vectr>> c1=r1/b1, c2=(r2-a2r1)/(b2-a2k1),c3=..etc
      
      g[i]=(c[i]-z[i-1]*g[i-1])/(x[i]-z[i-1]*f[i-1]) 
      
      #equating modified const vector w/ final soln vector
      
      d[v-1]=g[v-1]

#d stores the final solution for each elm simultaneously by backward subst

'Final Backward Substitution'

for i in range(v-1,0,-1):
      
      d[i-1]=g[i-1]-f[i-1]*d[i]      #backw subst (x1=r1-x2*k1, x2=r2-""..)
      
T=d                                  #soln by Thomas' algorithm          
W=np.linalg.solve(M,c)               #soln by in-built funcn                             

print("The solution vector by inbuilt numpy func: ",W)  

print("The solution vector by Thomas' Algorithm:",d)   

'PLOTTING'

#for plotting Y, we create a zero vectr w/ soln values in between

Y=np.zeros(n+1)
for i in range(1,n):
    Y[i]=d[i-1]

#for plotting X, we create a vector to store subsequent values in a loop
    
X=[]
for i in range(n+1):
        x=i*h
        X.append(x)
                 
ourplot.title("FINITE DIFFERENCE GRAPH")
ourplot.xlabel("X")
ourplot.ylabel("Plotting discrete y-values")
ourplot.plot(X,Y)
ourplot.show()