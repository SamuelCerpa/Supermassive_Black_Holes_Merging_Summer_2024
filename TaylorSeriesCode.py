import math
import numpy as np
#Taylor series of exponential and trigonometric functions

#Determine the factorial of any number 'n' then 'n!' if 'n!=n*(n-1)*(n-2)...*1, when 'n' is a natural number or zero
#Code on python
v1=input("Write the number of the factorial:") #Insert data
n=float(v1)
i=1
f=1 #f=1 allows '0!'
while(i<=n):
    f=f*i
    i=i+1
    # Write 'print(f)' to show the factorials for each number less than 'n' except '0!'
print("The factorial of",n,"is",f) #Show data

#Determine the value of any 'e^x' where 'e' is euler number and 'x' is any real number
#As of Taylor Series consider e^x= (x^0/0!)+(x^1/1!)+(x^2/2!)+....+(x^n/n!)
#Code on python
v2=input("Write the exponent of 'e':") #Insert data
x=float(v2) #'x' is the variable
n=1 # n is the power and factorial number increasing one by one 
i=1
f=1
w=1
while (n<=50): #The accuracy of 'e^x' is related to the inequality
    f=f*i
    i=i+1
    a=(x**n)/f #'a' is each expression of the serie '(x^n)/n!'
    w=w+a #'w' allows to sum each expression of 'a'
    n=n+1
print("'e' elevated to",x,"is",w) #Show data

#prototype
#Determine any 'sin(x)' where 'x' is on radians
#As of Taylor Series sin(x)=(x^1/1!)-(x^3/3!)+(x^5/5!)-...((-1)*(-1)^2n)*x^(2n-1)/(2n-1)!
v3=input("Write your angle:") #Write the angle on sexagesimal system
x=float(v3)*np.pi/180 #Convert for sexagesimal to radians
n=1
i=1
f=1
w=1
s=0
p=1
while(n<=8):
    f=f*i
    i=i+1
    a=(x**(n))
    if i % 2 == 0:
        w=(a/f)
    else:
        w=(0/f)
    n=n+1
    p=p+1
    print(a)
    print(f)
    print(w)
    