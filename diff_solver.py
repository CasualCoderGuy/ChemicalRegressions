#https://medium.com/@mathcube7/solving-differential-equations-analytically-with-python-b2dea629d50f
#https://docs.sympy.org/dev/guides/solving/solve-ode.html#use-the-solution-result


from sympy import *

init_printing()

AT = Function('AT')
ADT = Function('ADT')
ADM = Function('ADM')
AC = Function('AC')
t = symbols('t', real=True, positive=True)
k1, k2, k3, k4 = symbols('k1, k2, k3, k4', real=True, positive=True)

eq = [Eq(AT(t).diff(t), -k1 * AT(t)), Eq(ADT(t).diff(t), k1 * AT(t) - k2 * ADT(t)), Eq(ADM(t).diff(t), k2 * ADT(t) - k3 * ADM(t) + k4 * AC(t)), Eq(AC(t).diff(t), k3 * ADM(t) - k4 * AC(t))]
#difeq2 = Eq(ADT(t).diff(t), k1 * AT(t) - k2 * ADT(t))
#difeq2 = Eq(AD(t).diff(t), k1 * AT(t) - k2 * ADT(t))

AT0, ADT0, ADM0, AC0 = symbols('AT0, ADT0, ADM0, AC0')
initial = {
    AT0 : 100,
    ADT0 : 0,
    ADM0 : 0,
    AC0 : 0
}

sol = dsolve(eq, [AT(t), ADT(t), ADM(t), AC(t)])#.simplify()
print(sol)

initEq = Eq(sol[0].subs(t, 0).rhs.simplify(), AT0)
print(initEq)

consts = sol[0].atoms(Symbol).difference(eq[0].atoms(Symbol))
C1s = sorted(consts, key=lambda c: str(c))[0]
#coefs = solve([ini1, ini2], {C1, C2})
coefs = solve(initEq, C1s)
print(coefs)
C1 = symbols("C1")
coefsd = {C1: coefs}
print(sol[0].subs(coefsd))
