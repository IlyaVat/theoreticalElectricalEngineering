print("Hello world!")
import sys; print(sys.path)

import sys

import sympy
import mo

from sympy.parsing.sympy_parser import parse_expr
from sympy import laplace_transform

#sympy.pprint(mo.alap("1/(s**2+1)"))
sympy.pprint(mo.lap("10.0*Heaviside(t) "))
f=mo.sim("-11.4285714285714*hren((t-10.000000)) - 10.4481549985497*exp(-1.35355339059327*(t-10.000000))*hren((t-10.000000)) + 21.8767264271211*exp(-0.646446609406726*(t-10.000000))*hren((t-10.000000))")
f=mo.sep("(409600.0*s**2 + 4915200.0*s + 49152000.0)/(s*(128000.0*s**2 + 512000.0*s + 13312000.0))")

sympy.pprint(f);
print(str(sympy.nroots("(s*(128000.0*s**2 + 512000.0*s + 13312000.0))", n=30)))

#2.0/((4.0*s**2 + 8.0*s + 3.5))

f=parse_expr("0.5/((s**2 + 2.0*s + 0.875))")


x=sympy.Symbol("x",real=True)

f=f.subs({parse_expr("s"):x*parse_expr("I")})

print (f)
#print(mo.roo(f))
print("-----------------") 
print(mo.l_acx(f)) 
print("-----------------") 
print(mo.l_fcx(f))   

print(mo.alap("4.0*(s + 1)/(s*(2.0*s**2 + 13.0*s + 4.0))"))
"""
f=parse_expr("sin(t)")

sympy.pprint(sympy.simplify ('1'))
sympy.pprint(sympy.laplace_transform("sin(5*t)",'t','s'));

print("-----------------")

sympy.pprint(laplace_transform(f,sympy.Symbol('t'),sympy.Symbol('s')))
print(laplace_transform(f,sympy.Symbol('t'),sympy.Symbol('s'))[0])

x = sympy.Symbol('x')
y = sympy.Symbol('y')

r = sympy.Symbol ('r')

circle = sympy.Circle (sympy.Point (0, 0), r)
sympy.pprint("100")
sympy.pprint (circle.area)

"""