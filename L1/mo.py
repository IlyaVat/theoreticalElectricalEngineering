def hello():
    print('Hello, world!')

def my_measure(expr):
    from sympy import sqrt, simplify, count_ops, oo
    from sympy import Symbol, S
    
    count=0;
    # Discourage powers by giving POW a weight of 10
    count += count_ops(expr, visual=True).subs(Symbol('EXP'), 100)
    count += count_ops(expr, visual=True).subs(Symbol('HEAVISIDE'), 10)
    #count = count_ops(expr, visual=True).subs(Symbol('HEAVISIDE'), 1)
    
    # Every other operation gets a weight of 1 (the default)
    count = count.replace(Symbol, type(S.One))
    
    return count
def my_measure1(expr):
    from sympy import sqrt, simplify, count_ops, oo
    from sympy import Symbol, S
    count=0;
    # Discourage powers by giving POW a weight of 10
    #count += count_ops(expr, visual=True).subs(Symbol('HEAVISIDE'), 10)
    #count += count_ops(expr, visual=True).subs(Symbol('EXP'), 100)
    count += count_ops(expr, visual=True).subs(Symbol('NEG'), -1)
    
    # Every other operation gets a weight of 1 (the default)
    count = count.replace(Symbol, type(S.One))
    
    #print(count_ops(expr, visual=True))
    #print(count)
    
    return count

def sep(expr):
    import sympy
    return str(sympy.expand(expr))

def im(expr):
    import sympy
    return str(sympy.im(expr))
def re(expr):
    import sympy
    return str(sympy.re(expr))
def l_acx(expr):
    import sympy
    fi=sympy.im(expr)
    fr=sympy.re(expr)

    return str(sympy.simplify(sim(sympy.sqrt(fi*fi+fr*fr))))

def l_fcx(expr):
    import sympy
    fi=sympy.im(expr)
    fr=sympy.re(expr)

    return str(sympy.simplify(sympy.atan(sim(fi/fr))+sympy.Heaviside(-fr)*sympy.pi))

def roo(expr):
    import sympy
    from sympy import Poly
    print(expr)
    res=Poly(expr).all_roots()
    resok=[]
    for i in iter(res):
        resok.append(i.evalf())
    print(resok)
    return str(resok)

def sim(n):
    import sympy
    return str(sympy.simplify(n, measure=my_measure1, ratio=10))
def lap(n):
    import sympy
    return str(sympy.laplace_transform(n,'t','s')[0])
def alap(n):
    from sympy.integrals.transforms import inverse_laplace_transform
    from sympy import separatevars
    import sympy
    from sympy import simplify

    from sympy.parsing.sympy_parser import parse_expr

    x=sympy.Symbol("s",real=True)
    f=parse_expr(n)
    f=f.subs({parse_expr("s"):x})
    
    f=inverse_laplace_transform(f,x,'t')
    f=separatevars(f, force=True)
    #print(f)
    f=simplify(f, measure=my_measure1, ratio=10)
    
    return str(f)
def eva(n,n1,n2):
    import sympy
    from sympy.parsing.sympy_parser import parse_expr

    #sympy.pprint(parse_expr(n).xreplace({parse_expr(n1):parse_expr(n2)}))
    return str((parse_expr(n).xreplace({parse_expr(n1):parse_expr(n2)})).evalf())
def dif(n,n1):
    import sympy
    from sympy import sqrt, diff
    #from sympy.parsing.sympy_parser import parse_expr
    print(diff(n, n1,1))
    return str(diff(n, n1,1))
if __name__ == "__main__":
    hello()