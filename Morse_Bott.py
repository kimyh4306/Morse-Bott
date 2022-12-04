from itertools import chain, combinations, combinations_with_replacement, product
import numpy as np
from scipy import linalg
from operator import mul
from operator import itemgetter
from functools import reduce
from math import floor
from fractions import Fraction
import sympy
from sympy import poly
from sympy import Rational, nsimplify

def integ(f):
    return floor(f)

def frac(f):
    return Fraction(f,1)

def read_nr():
    temp = input("Enter coefficient:")
    try:
        number = int(temp)
    except ValueError:
        print("Integer coefficients")
        
    return number


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(0,len(s)+1))


def subsets(s):
    return list(map(set, powerset(s)) )

def tfconverter(w):
    temp=[]
    for i in range(0,len(w)):
      temp.append(1-w[i])
    return temp

def Singleorbit(a, m, n):
    I = []
    K = []
    temp = []
    for t in range(0,m):
      I.append(t)
    J = subsets(I)
    for t in J:
      for s in I:
        if s in t:
          temp.append(1)
        else:
          temp.append(0)
      K.append(temp)
      temp=[]
    #print(K)
    Single=[]
    tf = 0
    for t in K:
      for s in a:
        if np.matmul(np.array(s), np.array(tfconverter(t)).transpose()) == 0:
          tf=1
      if tf == 0:
        Single.append(t)  
      tf=0  
    del Single[0]
    return Single

def tempset(iterable, n):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range (n, n+1))

def rep_set(iterable, n):
    s = list(iterable)
    return chain.from_iterable(combinations_with_replacement(s, r) for r in range (n, n+1))


def length_subset(s, n):
    return list(map(set, tempset(s, n)) )

def rep_subset(s, n):
    return list(map(set, rep_set(s, n)) )


def gcd(a, b):
    if a < b:
        r = a
        a = b
        b = r
    r = b
    while( r != 0 ):
        r = a % b
        a = b
        b = r
    return a


def lcm(a, b):
    tmp = gcd(a, b)
    a = a // tmp
    return a * b


def gcd_set(a): #computes gcd of a set of integers
    nr_exp = len(a)
    if nr_exp == 0:
        return 1
    gcd_temp=a[0]
    for i in range(1, nr_exp):
        gcd_temp = gcd(a[i], gcd_temp)
    return gcd_temp


def lcm_set(a): #computes lcm of a set of integers
    nr_exp = len(a)
    prod = 1
    for i in range(0, nr_exp):
        prod = lcm(a[i], prod)
    return prod
    

def compare_value(a):
#Finds if all elements in a list are same.
#This function is also not used in the main code!
    if min(a) == max(a):
        return 1
    else:
        return 0
    
    
def sgn_odd_even(x):
    if x%2==0:
        return 1
    else:
        return -1
    

def transpose(a, m, n):
    if m > n:
        temp = [[0 for x in range(n)] for y in range(m)]
        for i in range (0,m):
            for j in range (0,n):
                temp[i][j]=a[j][i]
    else:
        temp = a
    return temp


def det_list(a, m, n):
    #Finds all nonzero minor determinants
    list=[]
    I=[]
    for t in range(n):
        I.append(t)
    ind=length_subset(I, m)
    for i in range(len(ind)):
        matrix_temp=[]
        for j in ind[i]:
            matrix_temp.append(a[j])
        mat=np.matrix(np.array(matrix_temp))
        if linalg.det(mat) != 0:
            list.append(abs(int(linalg.det(mat))))
    return list
 
def whp_detecter(a, m, n, d):
#Detects if coefficient rows lie on a hypersurface.
#Quick way to find if coefficients determine a weighted homogeneous polynomial (whp)
#a is the coefficient matrix
    if n == 1:
        return True
    else:
        temp = [[0 for x in range(m)] for y in range(n-1)]
        for i in range (0,n-1):
            for j in range (0,m):
                temp[i][j]=a[i+1][j] - a[0][j]
        ind = sympy.Matrix(np.array(temp)).nullspace()
        if len(ind) == 0:
            print("Polynomial is not weighted homogeneous")
            return False
        else:
            return True

"""        
def whp_finder(a, m, n, d): #Find weights by brute force repetition
#We do not use this function in the main code!
    I=[]
    for t in range(d):
        I.append(t+1)
    J=[]
    for t in range(n):
            J.append(t)
    large_set=[p for p in product(I, repeat=m)]
    for i in range(len(large_set)):
        for j in range(n):
            J[j]=np.inner(list(large_set[i]), a[j])
            if compare_value(J) == 1:
                print("Polynomial is weighted homogeneous of weights", large_set[i], " and total weight",J[0], ".")
                temper = list(large_set[i])
                temper.append(J[0])
                return temper
            else:
                if i==len(large_set)-1:
                    if j==n-1:
                        print("Polynomial is not weighted homogeneous.")
                        return []
"""
            
def whp_solver(Coeff, m, n, d):#Computes weights of whp
    B=[]
    for t in range(n):
        B.append(d)
    Bmat = np.matrix(np.array(B))#Converts Matrix to Numpy
    try:
        A = np.linalg.pinv(np.matrix(np.array(Coeff)))#Computes Pseudoinverse
    except np.linalg.LinAlgError:
        print("Polynomial is not weighted homogeneous.")
        return []
    x = np.matmul(A, Bmat.T)
    x = np.rint(x)
    y = x.tolist()
    weights = [0 for u in range(m)]
    for i in range (0,m):
        weights[i]=int(y[i][0])#Float->Int
    weights.append(d)#weights=(w_1,\cdots,w_n,d)
    gcd = gcd_set(weights)
    for i in range (0,m+1):
        weights[i]=weights[i] // gcd #Divide by gcd(w_1,\cdots,w_n,d)
    print("Polynomial is weighted homogeneous of weights", weights)
    return weights
                     
                        
def kappa(a,d):
    gcd_a=gcd_set(a)
    d = d // gcd_a
    for q in range(0,len(a)):
        a[q] = a[q] // gcd_a
    I=[]
    for t in range(0, len(a)):
        I.append(t)
    k=0
    sets=subsets(I)

    for t in range(0, len(sets)):
        a_s=[]
        for i in (sets[t]):
            a_s.append( a[i] )
        if t==0 :
            k = k + sgn_odd_even(len(I))
        else:
            prod_a_s = reduce(mul, a_s, 1)
            gcd_a_s = gcd_set(a_s)
            gcd_temp = gcd(gcd_a_s,d)
            sgn = sgn_odd_even(len(I)-len(a_s))
            change = sgn * (d ** (len(a_s)-1)) * Rational(gcd_temp , prod_a_s)
            k = k + change
    return k


def milnor_nr(a,d):
    b=[]
    for t in range (0,len(a)):
        b.append(Fraction(d, a[t]))
        #print(b[t])
    am = map(lambda x:x-1,b)
    return int(reduce(mul, am, 1))


def homology_fiber (a,d):
    Hom=[]
    Hom.append(1)
    for t in range (1,2*len(a)-1):
        Hom.append(0)
    
    Hom[len(a)-1]=milnor_nr(a,d)
    return Hom

def homology_rational (a,d):
    Hom=[]
    n=len(a)-1
    if n==0:
        Hom.append(1)
        Hom.append(1)
    elif n==1:
        Hom.append(int(kappa(a,d))+1)
        Hom.append(int(kappa(a,d))+1)
    else:
        Hom.append(1)
        for t in range (1,n-1):
            Hom.append(0)
        Hom.append(int(kappa(a,d)))
        Hom.append(int(kappa(a,d)))
        for t in range (n+2,2*n):
            Hom.append(0)
        Hom.append(1)
    
    
    return Hom

def homology_S1 (a,d):
    Hom=[]
    n=len(a)-1
    if n==1:
        Hom.append(int(kappa(a,d))+1)
    else:
        for t in range (0,2*n-1):
                if sgn_odd_even(t)==1:
                    Hom.append(1)
                else:
                    Hom.append(0)
        Hom[n-1]=Hom[n-1]+int(kappa(a,d))
    
    return Hom



def C_fct(a, Is, I ,d):#Randell's C Function
    if Is==[]:
        return Rational( d , lcm_set(a) )
    else:
        prod=1
        sets=subsets(Is)
        for It in sets:
            if len(It)!=len(Is):                
                prod = prod*C_fct(a, It, I, d)
        a_compl = []
        for i in I:
            if i not in Is:
                a_compl.append(a[i])
        return Rational ( Rational( d , lcm_set(a_compl) ) , prod )

    
def torsion(a, d):
    I = []
    for t in range(0, len(a)):
        I.append(t)
    d_list=[]
    #compute k-values and r
    k_vals=[]
    sets=subsets(I)
    r=0
    for Is in sets:
        #compute when odd, otherwise get 0
        if sgn_odd_even(len(I)-len(Is) ) == -1:
            a_s=[]
            for i in Is:
                a_s.append(a[i])            
            kappa_a_s=kappa(a_s, d)            
            k_vals.append(kappa_a_s)
            if kappa_a_s>r:
                r=int(kappa_a_s)
        else:
            k_vals.append(0)
    del k_vals[0]
    del sets[0]
    for j in range(1,r+1):
        prod=1
        counter=0
        for Is in sets:
            if k_vals[counter]>=j:
                prod=prod*C_fct(a, Is, I, d)
            counter=counter+1
        d_list.append(int(prod))
    return d_list


def matrix_trim(a, m, n):
#Detects same monomial appearing more than once
    a_mat = np.matrix(np.array(a))
    duplicate = np.where(~a_mat.any(axis=0))[1]
    print(duplicate)
    coeff_trim = [[0 for x in range(m - len(duplicate))] for y in range(n)] 
    temp_nr=0
    for j in range (0,m):
        if j in duplicate:
            temp_nr = temp_nr + 1
        else:
            for i in range (0,n):
                coeff_trim[i][j-temp_nr]=a[i][j]
    return(coeff_trim)

def ComputePeriodsAndOrbits(w, single):
    I = []
    K = []
    temp = []
    for t in range(0,len(w)):
        I.append(t)
    sets = subsets(I)
    for t in sets:
        for s in I:
            if s in t:
                temp.append(1)
            else:
                temp.append(0)
        K.append(temp)
        temp=[]
    period_list = []
    orbits = []
    per_list = []
        
    for t in range(0,len(sets)):        
        a_s=[]#a_s = weights
        for i in (sets[t]):
            a_s.append( w[i] )
            
        if K[t] in single:
            gcd_a_s=gcd_set(a_s)
            period_list.append( [Fraction(d , gcd_a_s), len(a_s), 1, a_s, sets[t]] )
                #put 1 as index at first
            orbits.append(a_s)
            per_list.append (Fraction(d, gcd_a_s))
            
        if len(a_s)>1:
            gcd_a_s=gcd_set(a_s)
            found=0
            for i in range(0, len(period_list) ):#finds largest MorseBott Component
                if (period_list[i])[0]== Fraction (d , gcd_a_s): 
                    found=1
                    if (period_list[i])[1] < len(a_s):
                        (period_list[i])[1] = len(a_s)
                        (period_list[i])[3] = a_s
                        (period_list[i])[4] = sets[t]
                        orbits[i] = a_s
            if found==0:
                period_list.append( [Fraction(d, gcd_a_s), len(a_s), 1, a_s, sets[t]] )
                #put 1 as index at first
                orbits.append(a_s)
                per_list.append (Fraction(d , gcd_a_s))
    return [I, a_s, period_list, orbits, per_list]

def MorseBottComponent(a, w, d, m, n):
    Comp = []
    denom = []
    
    periods = ComputePeriodsAndOrbits(w, Singleorbit(a, m, n))[2] #List of periods
    periods.sort(key=itemgetter(0)) #sort periods
    for s in range(0,len(periods)):
        denom.append(periods[s][0].denominator)
    mult = lcm_set(denom)  
    for i in range(d*mult):
        Comp.append([]) 
    I = []
    for t in range(0,len(w)):
        I.append(t)
    for i in range (len(periods)):
        for j in range (1, integ(Fraction(d, (periods[i][0]))) + 1):
            tem = None #refresh variable
            tem = j*(integ(frac(mult)*periods[i][0]))-1 #jth cover
            Comp[tem] = periods[i][:]#same as below code, but iteration problems
            #Comp[j*(periods[i][0])-1] = periods[i][:]
            if len(periods[i][3]) == 1:
                Comp[tem][2] = IndShift(w, d, I, list(periods[i][4]), periods[i][0], j)
            else:  
                Comp[tem][2] = IndShift(w, d, I, list(periods[i][4]), periods[i][0], j) - len(periods[i][3]) + 2
            #period_list[i][2] is always 1 initially
    for i in range (len(Comp)):
        if Comp[i] == []:
            continue
        else:
            a_s = Comp[i][3][:]
            a_s.append(d)
            component = a_s[:]
            gcd = gcd_set(component)
            for s in range(len(component)):
                component[s] = a_s[s] // gcd
            w = component[-1]
            del component[-1]
            Comp[i].append(homology_rational(component, w))
    return Comp

def MorseBottComponent_S1(a, w, d, m, n):
    Comp = []
    denom = []
    
    periods = ComputePeriodsAndOrbits(w, Singleorbit(a, m, n))[2] #List of periods
    periods.sort(key=itemgetter(0)) #sort periods
    for s in range(0,len(periods)):
        denom.append(periods[s][0].denominator)
    mult = lcm_set(denom)  
    for i in range(d*mult):
        Comp.append([]) 
    I = []
    for t in range(0,len(w)):
        I.append(t)
    for i in range (len(periods)):
        for j in range (1, integ(Fraction(d, (periods[i][0]))) + 1):
            tem = None #refresh variable
            tem = j*(integ(frac(mult)*periods[i][0]))-1 #jth cover
            Comp[tem] = periods[i][:]#same as below code, but iteration problems
            #Comp[j*(periods[i][0])-1] = periods[i][:]
            if len(periods[i][3]) == 1:
                Comp[tem][2] = IndShift(w, d, I, list(periods[i][4]), periods[i][0], j)
            else:  
                Comp[tem][2] = IndShift(w, d, I, list(periods[i][4]), periods[i][0], j) - len(periods[i][3]) + 2
            #period_list[i][2] is always 1 initially
    for i in range (len(Comp)):
        if Comp[i] == []:
            continue
        else:
            a_s = Comp[i][3][:]
            a_s.append(d)
            component = a_s[:]
            gcd = gcd_set(component)
            for s in range(len(component)):
                component[s] = a_s[s] // gcd
            w = component[-1]
            del component[-1]
            Comp[i].append(homology_S1(component, w))
    return Comp

def IndShift (w, d, I, Is, T, N):
    #Computes index shift
    temp = 0
    if (frac(2*N)*T).denominator == 1:
        temp = -integ(frac(2*N)*T)
    else:
        temp = -2*integ(frac(N)*T)-1

    for i in range(len(w)):
        if i in Is:
            temp = temp + integ(frac(2*N*w[i])*T) // d

        else:
            temp = temp + 2*(integ(frac(N*w[i])*T*Fraction(1,d)))  + 1
    return temp


def SSConverter (a):
#Rearranges Index of Spectral Sequence
    SS = []
    index = []
    temp = 0
    for i in range(len(a)):
        index.append(a[i][0] + a[i][1])
    mi = min(index)
    ma = max(index)
    for i in range(ma-mi+1):
        SS.append([])
    for i in range(len(a)):
        temp = a[i][0] + a[i][1]
        SS[temp-mi].append(a[i])
    return SS


def SS_rational (a):
#Converts Morse-Bott Components to E^1 page of Spectral Sequence
    m = len(a)
    SS = []
    p = 1
    for i in range(m):
        if a[i] == []:
            continue
        else:
            Hom = a[i][5]
            Shift = a[i][2]
            for j in range(len(Hom)):
                if Hom[j] == 0:
                    continue
                else:
                    SS.append([p,Shift-p+j,Hom[j]])
            p = p+1
    return SS

def Euler (a):
# Computes Mean Euler Characteristic
    chi = 0
    ind_max = 0
    ind_min = 0
    for i in range(len(a)):
        if len(a[i])==0:
            chi = chi
        else:
            if sgn_odd_even(a[i][0][0]+a[i][0][1])==1:
                for s in range(len(a[i])):
                    chi = chi + a[i][s][2]

            else:
                for s in range(len(a[i])):
                    chi = chi - a[i][s][2]
    ind_max=a[-1][0][0]+a[-1][0][1]
    ind_min=a[0][0][0]+a[0][0][1]
            
    return Fraction(chi, ind_max-ind_min)

print('Enter Polynomial')
m=int(input("Number of Variables: "))
n=int(input("Enter number of monomials in polynomial: "))
print("Enter coefficients for each monomial: ")

coeff = [[0 for x in range(m)] for y in range(n)] 
for i in range (0,n):
    for j in range (0,m):
        elt=read_nr()
        coeff[i][j]=elt
        print("------------------------------")
        if j == m-1:    
            print("Next Monomial")

    #coeff=[[3,0,0,0],[1,3,0,0],[0,0,2,0],[0,0,0,2]]
    #E7 Singularity for n+2=4 variables

coeff = matrix_trim(coeff,m,n)
m = len(coeff[0])
n = len(coeff)
print("Reduced coefficient matrix is: ", coeff)
    
coeff_temp = np.array(transpose(coeff, m, n))
dets=det_list(coeff_temp, min(m,n), max(m,n))
d_max=lcm_set(dets)
if whp_detecter(coeff,m,n,d_max) == True:
    weight = whp_solver(coeff,m,n,d_max)
    n = m
    d = weight[-1]
    del weight[-1]
    print("The Z-homology of the Milnor Fiber is:", homology_fiber(weight, d))
    print("The Q-homology of the link of singularity is:", homology_rational(weight,d))
    print("The middle dimension Z-torsion of the link of singularity is:",torsion(weight,d)) 
        #print("Period list: ")
        #for k in range(d):
            #print(MorseBottComponent(weight, d)[k])
    print()
    print("Spectral Sequence for positive homology has terms:")
    print()
    print(list(SSConverter(SS_rational(MorseBottComponent(coeff, weight, d, m, n)))))
    print()
    print("Mean Euler Characteristic for S^1 Equivariant SH is:")
    print()
    print(list(SSConverter(SS_rational(MorseBottComponent_S1(coeff, weight, d, m, n)))))
    print()
    print(Euler(list(SSConverter(SS_rational(MorseBottComponent_S1(coeff, weight, d, m, n))))))
