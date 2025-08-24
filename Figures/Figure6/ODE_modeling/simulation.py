import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import pandas as pd
from multiprocessing import Pool, cpu_count
from datetime import datetime
import time

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


def track_switching(P0, params, midpoint, signal_params):
    """Track switching dynamics with adaptive time steps"""
    signal_strength, signal_duration = signal_params
    params_with_signal = (*params[:6], signal_strength)
    
    # Three-phase time evaluation for accuracy
    t_eval = np.concatenate([
        np.linspace(0, 1, 50),  # Very dense early
        np.linspace(1.1, min(10, signal_duration)),  # Medium resolution
        np.linspace(min(10, signal_duration)+0.1, signal_duration)  # Coarse late
    ])

    try:
        P = odeint(model, P0, t_eval, args=(params_with_signal,)).flatten()
        crossing_idx = np.where(P > midpoint)[0]
        
        if len(crossing_idx) > 0:
            return t_eval[crossing_idx[0]], P[-1]  # Return time and final P
        return MAX_SWITCH_TIME, P[-1]  # Didn't cross midpoint
    except:
        return np.nan, np.nan  # Simulation failed


def simulate_cell_pair(shared_params, β0_low, β0_high, signal_params, t1, t2):
    α, β1, K_d, δ, n = shared_params
    signal_strength, signal_duration = signal_params
    
    # Find true initial stable states (without signal)
    params_no_signal_low = (α, β0_low, β1, K_d, δ, n, 0)
    params_no_signal_high = (α, β0_high, β1, K_d, δ, n, 0)
    
    stable_states_low = find_stable_states(params_no_signal_low)
    stable_states_high = find_stable_states(params_no_signal_high)

 #   print (stable_states_low)
  #  print (stable_states_high)
    
    # Use smallest stable state as initial condition
    P0_low = stable_states_low[0] if len(stable_states_low) == 3 else (α * β0_low / δ)
    P0_high = stable_states_high[0] if len(stable_states_high) == 3 else (α * β0_high / δ)
    
    # Phase 1: With signal
    params_signal_low = (α, β0_low, β1, K_d, δ, n, signal_strength)
    params_signal_high = (α, β0_high, β1, K_d, δ, n, signal_strength)
    
    t_phase1 = np.linspace(0, signal_duration, 100)
    P_mid_low = odeint(model, P0_low, t_phase1, args=(params_signal_low,))[-1, 0]
    P_mid_high = odeint(model, P0_high, t_phase1, args=(params_signal_high,))[-1, 0]
    
    # Phase 2: Remove signal
    t_phase2 = np.linspace(0, t2, 100)
    P_final_low = odeint(model, P_mid_low, t_phase2, args=(params_no_signal_low,))[-1, 0]
    P_final_high = odeint(model, P_mid_high, t_phase2, args=(params_no_signal_high,))[-1, 0]


    
    # Reference high state (approximate)
 #   high_state_ref = α*(β0_high + β1)/δ if δ > 0 else 100
    
    # Determine if switching occurred (using 75% of reference high state as threshold)
 #   threshold = 0.75 * high_state_ref
 #   switched_low = P_final_low > threshold
 #   switched_high = P_final_high > threshold
    switched_low = False
    stable_state_low_1 = 0
    stable_state_low_2 = 0
    stable_state_low_3 = 0
    switched_time_low = 0
    if len(stable_states_low) == 3:
        mid_point_low = stable_states_low[1]
        if ( P_mid_low > mid_point_low and P_final_low > mid_point_low ):
            switched_low = True
            switched_time_low = track_switching(P0_low, params_signal_low, mid_point_low, signal_params )[0]
        stable_state_low_1, stable_state_low_2, stable_state_low_3 = stable_states_low
    elif len(stable_states_low) == 1:
        stable_state_low_1 = stable_states_low[0]
        

    switched_high = False
    switched_time_high = 0
    stable_state_high_1 = 0
    stable_state_high_2 = 0
    stable_state_high_3 = 0
    if len(stable_states_high) == 3:
        mid_point_high = stable_states_high[1]
        if ( P_mid_high > mid_point_high and P_final_high > mid_point_high ):
            switched_high = True
            switched_time_high = track_switching(P0_high, params_signal_high, mid_point_high, signal_params )[0]
        stable_state_high_1, stable_state_high_2, stable_state_high_3 = stable_states_high
    elif len(stable_states_high) == 1:
        stable_state_high_1 = stable_states_high[0]
    
    return {
        'beta0_low': β0_low,
        'beta0_high': β0_high,
        'initial_P_low': P0_low,
        'initial_P_high': P0_high,
        'mid_P_low': P_mid_low,
        'mid_P_high': P_mid_high,
        'final_P_low': P_final_low,
        'final_P_high': P_final_high,
        'switched_low': switched_low,
        'switched_time_low': switched_time_low,
        'switched_high': switched_high,
        'switched_time_high': switched_time_high,
        'num_states_low': len(stable_states_low),
        'num_states_high': len(stable_states_high),
    #    'stable_states_low': stable_states_low,
    #    'stable_states_high': stable_states_high,
        'stable_state_1_low': stable_state_low_1,
        'stable_state_2_low': stable_state_low_2,
        'stable_state_3_low': stable_state_low_3,
        'stable_state_1_high': stable_state_high_1,
        'stable_state_2_high': stable_state_high_2,
        'stable_state_3_high': stable_state_high_3,
        'signal_strength': signal_strength,
        'signal_duration': signal_duration,
        'deta': δ,
        'beta1': β1,
        'K_d': K_d,
        'n': n
    }



def run_test():
    α = 1.0
    β1 = 2.0
    K_d = 70
    δ = 0.01
    n = 5
    signal_strength = 1.0
    signal_duration = 50.0
    β0_low = 0.15
    β0_high = 0.2
    
    t1= signal_duration
    t2 = 1000

    Result = simulate_cell_pair(
        shared_params = (α, β1, K_d, δ, n),
        β0_low = β0_low,
        β0_high = β0_high,
        signal_params = (signal_strength, signal_duration),
        t1 = signal_duration,
        t2 = t2
    )
    print ( "Results for the test:\n" )
    print ("beta0_low: %f\n" %Result['beta0_low'])
    print ("beta0_high: %f\n" %Result['beta0_high'])
    print ("initial_P_low: %f\n" %Result['initial_P_low'])
    print ("initial_P_high: %f\n" %Result['initial_P_high'])
    print ("mid_P_low: %f\n" %Result['mid_P_low'])
    print ("mid_P_high: %f\n" %Result['mid_P_high'])
    print ("final_P_low: %f\n" %Result['final_P_low'])
    print ("final_P_high: %f\n" %Result['final_P_high'])
    print ("switched_low: ", Result['switched_low'], '\n' )
    print ("switched_time_low: ", Result['switched_time_low'], '\n' )
    print ("switched_high: ", Result['switched_high'], '\n' )
    print ("switched_time_high: ", Result['switched_time_high'], '\n' )
    print ("num_states_low: %d\n" % Result['num_states_low'])
    print ("num_states_high: %d\n" % Result['num_states_high'])
    print ("stable_states_low: ", Result['stable_state_1_low'],Result['stable_state_2_low'], Result['stable_state_3_low'],'\n')
    print ("stable_states_high: ", Result['stable_state_1_high'],Result['stable_state_2_high'], Result['stable_state_3_high'], '\n')
    print ("signal_strength: %f\n" % Result['signal_strength'])
    print ("signal_duration: %f\n" % Result['signal_duration'])
    



def run_comparative_study(num_simulations=1000):
    np.random.seed(42)
    
    # Carefully chosen parameter ranges
    β1_range = (0.1, 10)  # Wider range including very weak feedback
    K_d_range = (10, 100.0)  # Full biological range
    δ_range = (0.001, 0.05)  # Very slow to fast degradation
    n_choices = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    # Signal parameters that can fail
    signal_strength_range = (0.1, 10.0)  # Weak signals that may not trigger switching
    signal_duration_range = (1, 50.0)  # Short pulses
    
    # Ensure distinct basal rates
    β0_low_range = (0.001, 0.05)
    β0_high_range = (0.1, 0.5)
    
    results = []
    for _ in range(num_simulations):
        # Sample parameters log-uniformly where appropriate
        β1 = 10**np.random.uniform(np.log10(β1_range[0]), np.log10(β1_range[1]))
        K_d = np.random.uniform(K_d_range[0], K_d_range[1])
        δ = 10**np.random.uniform(np.log10(δ_range[0]), np.log10(δ_range[1]))
        n = np.random.choice(n_choices)
        
        β0_low = np.random.uniform(*β0_low_range)
        β0_high = np.random.uniform(*β0_high_range)
        
        signal_strength = 10**np.random.uniform(np.log10(signal_strength_range[0]), np.log10(signal_strength_range[1]))
        signal_duration = np.random.uniform(*signal_duration_range)
        t2 = 1000
        
        result = simulate_cell_pair(
            shared_params = (1.0, β1, K_d, δ, n),
            β0_low = β0_low,
            β0_high = β0_high,
            signal_params = (signal_strength, signal_duration),
            t1 = signal_duration,
            t2 = t2
        )
        results.append(result)
    
    return pd.DataFrame(results)

def run_comparative_study_2(args):
    i, total = args[0], args[1]
    if i % 10000 == 0:
        print(f"Running simulation {i}/{total}")
    
    # Carefully chosen parameter ranges
    β1_range = (0.1, 10)  # Wider range including very weak feedback
    K_d_range = (10, 100.0)  # Full biological range
    δ_range = (0.001, 0.05)  # Very slow to fast degradation
    n_choices = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    # Signal parameters that can fail
    signal_strength_range = (0.1, 10.0)  # Weak signals that may not trigger switching
    signal_duration_range = (1, 50.0)  # Short pulses
    
    # Ensure distinct basal rates
    β0_low_range = (0.001, 0.1)
 #   β0_high_range = (0.1, 0.5)
    
    
    # Sample parameters log-uniformly where appropriate
    β1 = 10**np.random.uniform(np.log10(β1_range[0]), np.log10(β1_range[1]))
    K_d = np.random.uniform(K_d_range[0], K_d_range[1])
    δ = 10**np.random.uniform(np.log10(δ_range[0]), np.log10(δ_range[1]))
    n = np.random.choice(n_choices)
    
    β0_low = np.random.uniform(*β0_low_range)
    # fix ratio. Change ratio to simulate ratio = 2, 3, 4, 5, 6, 7, 8
    ratio = 9
    β0_high = β0_low * ratio
    
    signal_strength = 10**np.random.uniform(np.log10(signal_strength_range[0]), np.log10(signal_strength_range[1]))
    signal_duration = np.random.uniform(*signal_duration_range)
    t2 = 1000
    
    result = simulate_cell_pair(
        shared_params = (1.0, β1, K_d, δ, n),
        β0_low = β0_low,
        β0_high = β0_high,
        signal_params = (signal_strength, signal_duration),
        t1 = signal_duration,
        t2 = t2
    )
    
    return result

def run_large_simulation(num_simulations=100000):
    start_time = time.time()
    print(f"Starting {num_simulations:,} simulations...")
    
    with Pool(processes=max(1, cpu_count()-1)) as pool:
        results = [r for r in pool.map(
            run_comparative_study_2,
            [(i, num_simulations) for i in range(num_simulations)]
        ) if r is not None]
    
    if results:
        df = pd.DataFrame(results)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"switching_kinetics_{len(df)}cases_{timestamp}.csv"
        df.to_csv(filename, index=False)

    print(f"\nTotal runtime: {(time.time()-start_time)/60:.1f} minutes")
    return df if results else pd.DataFrame()

#if __name__ == "__main__":
#    df = run_test()
# Analysis
if __name__ == "__main__":
#    df = run_comparative_study(1000)
#    df.to_csv("bistable_switching_results.csv", index=False)
    df = run_large_simulation(num_simulations=100000)
    
    # Analyze by regime
#    df['regime'] = np.select(
#        [
 #           (df['num_states_low'] > 1) | (df['num_states_high'] > 1),
 #           (df['β1'] > 2) & (df['K_d'] < 50)
 #       ],
  #      ['bistable', 'strong_feedback'],
 #       default='monostable'
 #   )
    
 #   print("Switching probabilities:")
 #   print(df.groupby('regime')[['switched_low', 'switched_high']].mean())
  #  print("\nCases where initial states differ:")
 #   print(df[df['initial_P_low'] != df['initial_P_high']].shape[0]/len(df))