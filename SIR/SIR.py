import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Condizioni iniziali
N = 1000
I0 = 1
R0 = 0
S0 = N - I0 - R0
y0 =(S0, I0, R0)

# Parametri del modello (tempo finale, beta e gamma)
tf = 200
beta = 0.2
gamma = 0.1

# Il modello SIR
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt


# Risoluzione solo nei punti interi
t = np.linspace(0, tf, tf+1)
sol = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = sol.T

# Visualizzazione
fig = plt.figure()
plt.plot(t, S, 'b', lw=2, label='S')
plt.plot(t, I, 'r', lw=2, label='I')
plt.plot(t, R, 'g', lw=2, label='R')

plt.legend()
