import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import ndimage


# 定义WS网络参数
N = 1000  # 节点数量
K = 4  # 平均度
p = 0.5  # 重连概率

# 初始化WS网络
ws_graph = nx.watts_strogatz_graph(N, K, p)

# 模拟传播过程
total_time_steps = 1000
intera = 1

beta0_range = np.linspace(0, 0.5, num=10)  # beta0范围
#np.round(beta0_range, 2)
delta_range = np.linspace(0, 9, num=10)  # delta范围
print(beta0_range)
heatmap_data = np.zeros((len(delta_range), len(beta0_range)))  # 用于存储感染人数比例的二维数组
for beta0_index, beta0 in enumerate(beta0_range):
    print(beta0)
    print(str(beta0_index / 10 * 100) + '%')
    for delta_index, delta in enumerate(delta_range):
        print(delta)
        TT = []  # 针对每个beta值，收集五次迭代后的数组
        for _ in range(intera):
            # 其他代码省略...
            node_state = {}
            # 0: S; 1: I; 2: R.
            for node in ws_graph.nodes:
                node_state[node] = 0
            node_state[random.choice(list(ws_graph.nodes))] = 1

            gamma = 0.4  # I->R
            alpha = 0.5  # R->S

            susceptible_count_t = []
            infected_count_t = []
            recovered_count_t = []
            susceptible_count = sum(value == 0 for value in node_state.values())
            infected_count = sum(value == 1 for value in node_state.values())
            recovered_count = sum(value == 2 for value in node_state.values())
            # 计算初始态sir个数
            susceptible_count_t.append(susceptible_count)
            infected_count_t.append(infected_count)
            recovered_count_t.append(recovered_count)

            for t in range(total_time_steps):
                # 其他代码省略...
                for node in ws_graph.nodes:
                    if node_state[node] == 0:
                        # 动态感染率
                        sigma = 0.001  # 扩散率
                        miu = 0  # 漂移率
                        # 布朗运动模拟
                        dt = 1  # 时间步长
                        # 绘制传染率随时间变化的图像
                        dB = np.sqrt(dt) * np.random.normal(size=total_time_steps)  # 取size=total_time_steps个数，Brownian motion increment,按正态概率正负无穷取值
                        B = np.cumsum(dB)  # Brownian motion path，累加过后的数组
                        beta = beta0 * np.exp((miu - 0.5 * sigma ** 2) * t + sigma * B)  # beta是一个数组，为B是一个数组，SIR模型中的传染率β的随机漂移，
                        # 自我隔离强度
                        neis = list(ws_graph.neighbors(node))
                        total_xy_ratio = 1
                        for neighbor1 in neis:  # 遍历一代邻居
                            if node_state[neighbor1] == 1:
                                neigh_num = 0
                                infected_neis_count = 0
                                for neighbor2 in ws_graph.neighbors(neighbor1):  # 遍历一代邻居的邻居
                                    neigh_num += 1
                                    if node_state[neighbor2] == 1:
                                        infected_neis_count += 1
                                xy_ratio = 1 - (np.exp(- delta * infected_neis_count / neigh_num) * beta[t])
                                total_xy_ratio *= xy_ratio  # 将每个感染邻居的（1-cβ）值相乘
                            else:
                                continue
                        probai = 1 - total_xy_ratio

                        if 0 < random.uniform(0, 1) < probai:
                            node_state[node] = 1

                    elif node_state[node] == 1:
                        if 0 < random.uniform(0, 1) < gamma:
                            node_state[node] = 2

                    elif node_state[node] == 2:
                        if 0 < random.uniform(0, 1) < alpha:
                            node_state[node] = 0

                # print(f"计算{t}时刻sir态个体数目并保存在数组中")
                susceptible_count = sum(value == 0 for value in node_state.values())
                infected_count = sum(value == 1 for value in node_state.values())
                recovered_count = sum(value == 2 for value in node_state.values())

                susceptible_count_t.append(susceptible_count)
                infected_count_t.append(infected_count)
                recovered_count_t.append(recovered_count)

                f_infected_count_t = [value / N for value in infected_count_t]
                f_susceptible_count_t = [value / N for value in susceptible_count_t]
                f_recovered_count_t = [value / N for value in recovered_count_t]

            TT.append(f_infected_count_t)

        result_array = [sum(elements) / intera for elements in zip(*TT)]  # 将重复试验得到的数组的元素取均值，得到新数组
        heatmap_data[delta_index, beta0_index] = result_array[-1]  # 存储感染人数比例，-1表示最后一个时间步的数据

# 保存 heatmap_data 到一个 NumPy 文件
np.save('SIRS heatmap_data.npy', heatmap_data)
plt.imshow(heatmap_data, cmap='jet', aspect='equal', origin='lower', interpolation='gaussian')
#plt.title('Smoothed Heatmap')
plt.xlabel('β(0)')
plt.ylabel('δ')
plt.colorbar()

#高斯平滑
smoothed_heatmap_data = ndimage.uniform_filter(heatmap_data, size=4)

plt.subplot(1, 2, 1)
plt.imshow(smoothed_heatmap_data, cmap='jet', aspect='auto', origin='lower')
plt.title('Smoothed Heatmap (Gaussian)')
#plt.xticks(np.linspace(0, 0.5, len(a)), a)  # 使用自定义刻度a

plt.colorbar(label='I1')
plt.xlabel('beta0')
plt.ylabel('delta')


# 绘制原始热力图
plt.subplot(1, 2, 2)
plt.imshow(heatmap_data, cmap='jet', aspect='auto', origin='lower')
plt.title('Original Heatmap')
# 添加颜色条
plt.colorbar(label='I2')
# 添加坐标轴标签
plt.xlabel('beta0')
plt.ylabel('delta')
#plt.savefig('relit.pdf', format='pdf')

# 显示图形
plt.tight_layout()
plt.show()