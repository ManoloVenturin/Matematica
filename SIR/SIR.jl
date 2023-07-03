using DifferentialEquations
using Plots

# Condizioni iniziali
N = 1000
I0 = 1
R0 = 0
S0 = N - I0 - R0

# Parametri del modello (tempo finale, beta e gamma)
tf = 200
βvalue = 0.2
γvalue = 0.1

# Il modello SIR
function SIR(y, p, t)
    S, I, R = y
    β, γ = p
    dSdt = -β * S * I / N
    dIdt = β * S * I / N - γ * I
    dRdt = γ * I
    return [dSdt, dIdt, dRdt]
end

# Solution
y0 = [S0, I0, R0]
tspan = (0, tf)
params = [βvalue, γvalue]

prob = ODEProblem(SIR, y0, tspan, params)
sol = solve(prob)

# Plot
plot(sol)
# plot(sol.t, sol[1,:], label='S')
# plot!(sol.t, sol[2,:], label='I')
# plot!(sol.t, sol[3,:], label='R')



# Versione con il linguaggio di modellazione
using ModelingToolkit

# Definizione delle variabili
@parameters t β γ

# Variabili
@variables S(t) I(t) R(t)

# Operatore di derivazione temporale
Dt = Differential(t)

# Sistema di equaioni
eqs = [Dt(S) ~ -β * S * I / N,
       Dt(I) ~ β * S * I / N - γ * I,
       Dt(R) ~ γ * I]

@named de = ODESystem(eqs)

SIRmodel = ODEFunction(de, [S,I,R], [β,γ])

# Solution
y0 = [S0, I0, R0]
tspan = (0, tf)
params = [βvalue, γvalue]

prob2 = ODEProblem(SIRmodel, y0, tspan, params)
sol2 = solve(prob2)
plot(sol2)




