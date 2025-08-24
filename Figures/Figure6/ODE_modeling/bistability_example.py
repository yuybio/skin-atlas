import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def model(P, t, params):
    α, β0, β1, K_d, δ, n, signal = params
    activation = (P**n) / (K_d**n + P**n)
    dPdt = α * (β0 + β1*activation + signal) - δ*P
    return dPdt

def find_stable_states(params, num_guesses=100):
    """Find all stable states using multiple initial guesses"""
    def equation(P):
        return model(P, 0, params)  # dt=0 for steady state
    
#    print ("params: ", params)
    solutions = []
    guesses = np.linspace(0, 1000, num_guesses)  # From 1 to 1000
    tolerance = 1
    
    for guess in guesses:
        
        sol = fsolve(equation, guess)[0]
        if sol >=0 and equation(sol) < 1e-4 and equation(sol) > -1e-4:

            if not any(abs(sol - existing) < tolerance for existing in solutions):
                solutions.append(sol)
           
    return sorted(solutions)

# Parameters

α = 1.0
β0_low = 0.1
β0_high = 0.4
β1 = 2.9359349632193927
δ = 0.047755403396725495
K_d = 35.66448806808264
n = 6
signal_duration = 30
signal_strength = 1.2521681107514169
t2 = 100

params_no_signal_low = (α, β0_low, β1, K_d, δ, n, 0)
params_no_signal_high = (α, β0_high, β1, K_d, δ, n, 0)
    
stable_states_low = find_stable_states(params_no_signal_low)
stable_states_high = find_stable_states(params_no_signal_high)
print ("stable states:", stable_states_low, stable_states_high)

P0_low = stable_states_low[0] if len(stable_states_low) == 3 else (α * β0_low / δ)
P0_high = stable_states_high[0] if len(stable_states_high) == 3 else (α * β0_high / δ)
print (P0_low, P0_high)

# Phase 1: With signal
params_signal_low = (α, β0_low, β1, K_d, δ, n, signal_strength)
params_signal_high = (α, β0_high, β1, K_d, δ, n, signal_strength)
    
t_phase1 = np.linspace(0, signal_duration, 100)
P_phase1_low = odeint(model, P0_low, t_phase1, args=(params_signal_low,))
P_mid_low = P_phase1_low[-1, 0]
P_phase1_high = odeint(model, P0_high, t_phase1, args=(params_signal_high,))
P_mid_high = P_phase1_high[-1, 0]
print (P_mid_low, P_mid_high)
    
    # Phase 2: Remove signal
t_phase2 = np.linspace(0, t2, 100)
P_phase2_low = odeint(model, P_mid_low, t_phase2, args=(params_no_signal_low,))
P_final_low = P_phase2_low[-1,0]
P_phase2_high = odeint(model, P_mid_high, t_phase2, args=(params_no_signal_high,))
P_final_high = P_phase2_high[-1, 0]
print( P_final_low, P_final_high)

plt.figure(figsize=(7.5, 4.5))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 'black']
plt.plot(t_phase1, P_phase1_low, 'b', c= colors[0], label=f'X_low phase1')
plt.plot((t_phase2 + signal_duration), P_phase2_low, 'b', c= colors[1], label=f'X_low phase2')
plt.plot(t_phase1 , P_phase1_high, 'r', c= colors[2], label=f'X_high phase1')
plt.plot((t_phase2 + signal_duration), P_phase2_high, 'r', c= colors[3],  label=f'X_high phase2')

plt.xlabel('Time (s)'); plt.ylabel('P (nM)')
plt.title('Bistability in a Self-Activating Gene Circuit')
plt.legend(); plt.grid()
#plt.show()
plt.savefig('an_example.pdf')
