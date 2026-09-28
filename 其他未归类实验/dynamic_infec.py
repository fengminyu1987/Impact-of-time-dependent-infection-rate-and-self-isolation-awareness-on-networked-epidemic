import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
intera = 1000
total_time_steps = 50
# 动态感染率
sigma = 0.02  # 扩散率
miu = -0.01  # 漂移率
# 布朗运动模拟
dt = 1  # 时间步长

TT=[]

beta = [0] * (total_time_steps+1)

beta0_array = [0.1, 0.3, 0.5, 0.7, 0.9, 0.9,0.9,0.9,0.9,0.9]
final_result_array=[]
for beta0 in beta0_array:
    print(beta0)
    TT=[]
    for i in range(intera):
        beta[0]=beta0
        # 绘制传染率随时间变化的图像
        dB = np.sqrt(dt) * np.random.normal(
            size=total_time_steps)  # 取size=total_time_steps个数，Brownian motion increment,按正态概率正负无穷取值
        B = np.cumsum(dB)  # Brownian motion path，累加过后的数组
        for t in range(total_time_steps):
          beta[t+1] = (beta0 * np.exp((miu - 0.5 * sigma ** 2) * t + sigma * B[t]))  # SIR模型中的传染率β的随机漂移

        TT.append(beta)
    result_array = [sum(elements) / intera for elements in zip(*TT)]  # 将重复试验得到的数组的元素取均值，得到新数组
    final_result_array.append(result_array)  # 元素个数等于beta数组个数

for i in range(len(final_result_array)):
    plt.plot(range(total_time_steps+1), final_result_array[i], label="δ={}".format(beta0_array[i]))



plt.xlabel('Time')
plt.ylabel('BT')
#plt.xscale('log')
plt.xlim(0, 50)
plt.yticks(np.arange(0, 1.51, 0.1))
plt.legend() # 添加图例，即用来解释曲线代表的东西
plt.title(f'Time Evolution')
plt.grid()  # 添加网格线

plt.show()