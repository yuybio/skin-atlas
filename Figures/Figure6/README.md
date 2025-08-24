# Figure 6 — ODE Modeling

This folder contains python scripts used to **run ODE model simulation** for the Figure 6 panels.

---

## Contents

- `simulation.py`
- `stable_states_feedback.py`
- `stable_state_kbasal.py` 
- `bistability_example.py` 

## Dependencies

- **Python ≥ 3.9**
- `numpy`, `scipy` (ODE solvers via `scipy.integrate.solve_ivp`), `matplotlib`
- Optional: `pandas`, `seaborn`

---

## Quickstart

### 1) Run simulations on 100,000 randomized parameters with fixed ratio between X_high and X_low:

```bash
python simulation.py
```

> Edit Line 278 of the script to change ratios.

### 2) An example showing steady states VS. feedback strength.

```bash
Python steady_state_feedback.py
```

### 3) An example showing steady states VS. basal level.

```bash
Python steady_state_kbasal.py
```

### 4) An example showing bi-stability

```bash
Python bistability_example.py
```
