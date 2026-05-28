import sympy as sp
from sympy import Symbol, Function

# EDO 2O na forma:
# y'' + Py' + Qy = G        (01)
# Pode ser resolvida pelo MCFI em caso de P e Q constantes

# Edo Fácil:

# y'' + y = 3sen(2t) + tcos(2t)
# a = -5/9
# b = -1/3

# y'' + 2y' + 2y = 3e^-t + 2e^(-t)cos(t) + 4e^(-t)t²sen(t)

# Definição dos símbolos na equação original
t = Symbol('t')

a = Symbol('a')
b = Symbol('b')

# Definição das funções que compões G(t)
f = sp.sin(2*t)
g = sp.cos(2*t)*t

# Definição de P, Q e G na EDO de forma (01)
P = 0
Q = 1
G = 3*f + g

# Definição de y_ (y suposto)
y_ = a*f + b*g

# Montagem da EDO
edo = sp.Eq( sp.diff(y_, t, 2) + P*sp.diff(y_, t, 1) + Q*y_ , G)

# Resolução da EDO (obtendo sistema de equação)
edo_expandida = sp.expand(edo.lhs - edo.rhs)
termos_base = (f, g)

coeficientes = sp.collect(edo_expandida, termos_base, evaluate=False)
sistema = [sp.Eq(coef, 0) for coef in coeficientes.values()]

solucao = sp.linsolve(sistema, (a, b))
valores = list(solucao)[0]

va = valores[0]
vb = valores[1]

y = y_.subs({a:va, b:vb})
print(y)