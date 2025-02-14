from sympy import dsolve, symbols, Function, Eq, sin

x= symbols('x')

f=Function('f')

t=Eq(f(x).diff(x,1)-f(x)- sin(x),0)

print(dsolve(t,f(x)))
