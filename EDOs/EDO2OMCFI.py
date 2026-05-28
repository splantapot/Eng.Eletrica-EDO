import sympy as sp
from sympy import Symbol, Function

# EDO 2O na forma:
# y'' + Py' + Qy = G        (01)
# Pode ser resolvida pelo MCFI em caso de P e Q constantes

# y'' + 2y' + 2y = 3e^-t + 2e^(-t)cos(t) + 4e^(-t)t²sen(t)      (02)

# Definição dos símbolos na equação original
t = Symbol('t')

funcoes = [
    # Primeiro fator
    sp.exp(-t),
    # Terceiro fator
    (t**3)*sp.exp(-t)*sp.sin(t),
    (t**3)*sp.exp(-t)*sp.cos(t),
    (t**2)*sp.exp(-t)*sp.sin(t),
    (t**2)*sp.exp(-t)*sp.cos(t),
    (t)*sp.exp(-t)*sp.sin(t),
    (t)*sp.exp(-t)*sp.cos(t),
    # Segundo fator
    sp.exp(-t)*sp.sin(t),
    sp.exp(-t)*sp.cos(t)
]

coef_e_funcoes = [[Symbol(f"A{i}"), f] for i, f in enumerate(funcoes)]

# Definição de P, Q e G na EDO de forma (01)
P = 2
Q = 2
G = 2*coef_e_funcoes[0][1] + 2*coef_e_funcoes[8][1] + coef_e_funcoes[3][1]

# Definição de y_ (y suposto)
y_ = 0
for symbol, func in coef_e_funcoes:
    y_ = y_ + symbol*func

# Montagem da EDO
edo = sp.Eq( sp.diff(y_, t, 2) + P*sp.diff(y_, t, 1) + Q*y_ , G)

# Resolução da EDO (obtendo sistema de equação)
edo_expandida = sp.expand(edo.lhs - edo.rhs)
termos_base = tuple(f for s, f in coef_e_funcoes)

coeficientes = sp.collect(edo_expandida, termos_base, evaluate=False)
sistema = [sp.Eq(coef, 0) for coef in coeficientes.values()]

if (sistema.count(False)):
    print("Erro ao resolver a EDO. O Sistema possui uma equaçao falsa.")
else:
    solucao = sp.linsolve(sistema, tuple(s for s, f in coef_e_funcoes))
    # print(solucao)
    valores = list(solucao)[0]
    print("Coeficientes An:\n")
    print(valores)

    # x = Symbol('x')
    substituicoes = {} #{t:x}
    for coef, valor in zip(tuple(s for s, f in coef_e_funcoes), valores):
        if (valor == coef):
            substituicoes[coef] = 0
        else:
            substituicoes[coef] = valor

    print("\nSolução:\n")
    y = y_.subs(substituicoes)
    print(y)