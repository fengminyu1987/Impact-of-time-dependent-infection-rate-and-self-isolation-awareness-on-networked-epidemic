import numpy as np
import matplotlib.pyplot as plt

total_time_steps = 1000
sigma = 0.02
miu = -0.005
dt = 1

dB = np.sqrt(dt) * np.random.normal(size=total_time_steps)
B = np.cumsum(dB)
beta = np.zeros(total_time_steps)

beta0_array = [0.5, 0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]
for beta0 in beta0_array:
    for t in range(total_time_steps):

        beta[t] = beta0 * np.exp((miu - 0.5 * sigma ** 2) * t + sigma * B[t])
        print(beta[0])

    plt.plot(beta, label=f'beta0={beta0}')

plt.xlabel('Time')
plt.ylabel('BT')
#plt.xscale('log')
plt.xlim(0, 500)
plt.yticks(np.arange(0, 1.51, 0.1))
plt.legend() # 添加图例，即用来解释曲线代表的东西
plt.title(f'Time Evolution')
plt.grid()  # 添加网格线


plt.show()