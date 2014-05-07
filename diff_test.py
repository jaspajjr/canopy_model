import time
import numpy as np 
from sympy import *

start = time.time()
c, b, m = symbols("c, b, m")
z, n, p, q = symbols("z, n, p, q")

x = symbols("x")

g = (c / (1 + exp(-b * (x - m)))) * (z + (c * exp(-exp(-p * (x - n)))))

f = (0.8 / (1 + exp(-0.001 * (x - 2000)))) * (0.16 + (-0.8 * exp(-exp(-0.001 * (x - 700)))))

print type(f)
#f = exp(x/(1 + x))
fx = diff(f, x)
gx = integrate(g, x)
Fx = lambdify(x, fx, 'numpy')
Gx = lambdify(x, fx, 'numpy')
xi = np.linspace(0, 2, 5)
print Fx(xi)
print Gx(xi)

elapsed = time.time() - start
print elapsed