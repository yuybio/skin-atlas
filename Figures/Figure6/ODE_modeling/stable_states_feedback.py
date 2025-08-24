import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 定义稳态方程（dP/dt = 0）
def model(P, t, params):
    α, β0, β1, K_d, δ, n = params
    dPdt = α * (β0 + β1 * (P**n) / (K_d**n + P**n)) - δ * P
    return dPdt


# Parameters
#params = [1.0, 0.1, 3, 100, 0.01, 9]  # α, β0, β1, K_d, δ, n
t = np.linspace(0, 2000, 1000)         # Time (s)

# Define the model equation to solve (dP/dt = 0)
def steady_state_eq(P, params):
    α, β0, β1, K_d, δ, n = params
    return α * (β0 + β1 * (P**n) / (K_d**n + P**n)) - δ * P

# Find roots numerically
def find_steady_states(params, P_range):
    # Create a grid of initial guesses
    P_guesses = np.linspace(P_range[0], P_range[-1], 100)
    
    # Find unique roots
    solutions = []
    tolerance = 1
    
    for guess in P_guesses:
        sol = fsolve(steady_state_eq, guess, args=(params,))
        if steady_state_eq(sol[0], params) > 1e-4 or steady_state_eq(sol[0], params) < -1e-4:
            continue
        if sol[0] >= 0:  # Physical solutions only
            
            # Check if solution is new
            if not any(abs(sol[0] - existing) < tolerance for existing in solutions):
                solutions.append(sol[0])
    
    return np.sort(np.unique(np.round(solutions, 4)))

# 参数设置
n = 2       # Hill系数
K_d = 100.0      # 半饱和常数
delta = 0.05   # 降解速率
alpha = 1.0   # 生产系数

# 分岔分析参数
k_basal_list = [0.5]  # 不同基础激活水平
k_feedback_range = np.linspace(0, 40, 1000)  # 反馈强度扫描范围
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 'black']

plt.figure(figsize=(7.5, 5.5))


P_range = np.linspace(0, 1000, 100)

# 遍历不同k_basal值
for idx, k_basal in enumerate(k_basal_list):
    stable_high = []
    stable_middle = []
    stable_low = []
    
    # 扫描反馈强度参数
    for k_feedback in k_feedback_range:
        # 使用不同初始猜测寻找多解

        params = (alpha, k_basal, k_feedback, K_d, delta, n)
        solutions = find_steady_states(params, P_range)
        
        if ( len(solutions) == 3 ):
            stable_low.append( (k_feedback, solutions[0]))
            stable_middle.append( (k_feedback, solutions[1]))
            stable_high.append( (k_feedback, solutions[2]))
        elif ( len(solutions ) == 1 ):
            if ( solutions[0] < 100 ):
                stable_low.append( (k_feedback, solutions[0]))
            elif ( solutions[0] > 100 ):
                stable_high.append((k_feedback, solutions[0]) )
        
      
    
    # 绘制分岔曲线
    if stable_low:
        x, y = zip(*stable_low)
    #    plt.plot(x, y, color=colors[0], label=f'k_basal={k_basal}', linewidth=2)
        plt.plot(x, y, color=colors[0])
    if stable_middle:
        x, y = zip(*stable_middle)
        plt.plot(x, y, color=colors[1], linestyle='--')
#        plt.plot(x, y, color=colors[idx], linestyle='--' if k_basal>0.1 else '-')
    if stable_high:
        x, y = zip(*stable_high)
        plt.plot(x, y, color=colors[2] )
 #       plt.plot(x, y, color=colors[idx], label=f'k_basal={k_basal}', linewidth=2)

# 图形标注
plt.xlabel('Feedback Strength (k_feedback)', fontsize=12)
plt.ylabel('Steady State Concentration (P)', fontsize=12)
plt.title('Bifurcation Diagram of Positive Feedback System', fontsize=14)
plt.grid(alpha=0.3)
#plt.legend()
plt.xlim(5, 17.5)
plt.ylim(0, 350)

# 标注阈值区域
#plt.fill_between([2.8, 5], [0, 0], [6, 6], color='red', alpha=0.1)
#plt.text(3.5, 1.5, 'Bistable Region\n(Tumor possible)', ha='center', color='darkred')

#plt.show()
#plt.savefig('steady_state_feedback.svg')
plt.savefig('steady_state_feedback.pdf')
