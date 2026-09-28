#import Gaussian as Gaussian
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import ndimage
from matplotlib.ticker import LinearLocator

# 定义WS网络参数
N = 1000  # 节点数量
K = 4  # 平均度
p = 0.5  # 重连概率

# 初始化WS网络
ws_graph = nx.watts_strogatz_graph(N, K, p)

# 模拟传播过程
total_time_steps = 1000
intera = 10
#beta0_range = [0, 5/95, 10/90, 15/85, 20/80, 25/75, 30/70, 35/65, 40/60, 45/55, 50/50, 55/45, 60/40, 65/35, 70/30, 75/25, 80/20, 85/15, 90/10, 95/5, 100]  # beta0范围,20geshu
#beta0_range = [0, 10/90, 20/80, 30/70, 40/60, 50/50, 60/40, 70/30, 80/20, 90/10, 100]  # beta0范围,20geshu
beta0_range = np.linspace(0, 1, num = 22)
beta0_range = np.round(beta0_range, decimals=3)
print(beta0_range)

delta_range = np.linspace(0, 1, num = 22)  # delta范围,20geshu
delta_range = np.round(delta_range, decimals=3)
print("delta_range:", delta_range)

heatmap_data = np.zeros((len(beta0_range), len(delta_range)))  # 用于存储感染人数比例的二维数组
for delta_index, delta in enumerate(delta_range):
    print("delta:", delta)
    print(str(delta_index / 22 * 100) + '%')
    for beta0_index, beta0 in enumerate(beta0_range):
        print("beta0:", beta0)
        TT = []  # 针对每个beta值，收集五次迭代后的数组
        for _ in range(intera):
            # 其他代码省略...
            node_state = {}
            # 0: S; 1: I; 2: R.
            for node in ws_graph.nodes:
                node_state[node] = 0

            initial_infected_nodes = random.sample(list(ws_graph.nodes), 10)  # 随机选择10个节点
            for node in initial_infected_nodes:
                node_state[node] = 1

            gamma = 0.4  # I->R

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

            # 动态感染率
            sigma = 0.005  # 扩散率
            miu = -0.001  # 漂移率
            # 布朗运动模拟
            dt = 1  # 时间步长
            # 绘制传染率随时间变化的图像
            dB = np.sqrt(dt) * np.random.normal(size=total_time_steps)  # 取size=total_time_steps个数，Brownian motion increment,按正态概率正负无穷取值
            B = np.cumsum(dB)  # Brownian motion path，累加过后的数组
            beta = [0]*total_time_steps
            #每次重复试验，给出新的beta[t]  数组
            for t in range(total_time_steps):
                beta[t] = (beta0 * np.exp((miu - 0.5 * sigma ** 2) * t + sigma * B[t]))  # beta是一个数组，为B是一个数组，SIR模型中的传染率β的随机漂移，

                #试验进行1000步
            for t in range(total_time_steps):
                # 其他代码省略...
                for node in ws_graph.nodes:
                    #判断当前t每个节点的状态
                    if node_state[node] == 0:
                        # 自我隔离强度
                        neis = list(ws_graph.neighbors(node))
                        total_xy_ratio = 1
                        for neighbor1 in neis:  # 遍历一代邻居
                            if node_state[neighbor1] == 1:#若一代邻居不正常，计数 2代正常和感染人数，1-（1-beta）^n
                                neis2 = list(ws_graph.neighbors(neighbor1))
                                neigh_num = len(neis2) #2代邻居数
                                infected_neis_count = 0
                                for neighbor2 in neis2:  # 遍历2代邻居，找不正常数
                                    if node_state[neighbor2] == 1:
                                        infected_neis_count += 1#2代感染数
                                    else:
                                        continue
                                xy_ratio = 1 - (np.exp(- delta*infected_neis_count / neigh_num) * beta[t]) #没被2代感染概率
                                total_xy_ratio *= xy_ratio  # 总的没被感染率，每（1-cβ）值相乘
                                #print("打印 beta[t]")
                                #print(beta[t])
                            else:
                                continue#若某一代邻居正常，则不用管，判断下个一代邻居
                        probai = 1 - total_xy_ratio #若一代邻居都正常，probai=1-1=0

                        if 0 < random.uniform(0, 1) < probai:
                            node_state[node] = 1

                    elif node_state[node] == 1:
                        if 0 < random.uniform(0, 1) < gamma:
                            node_state[node] = 2

                # print(f"计算{t}时刻sir态个体数目并保存在数组中")
                susceptible_count = sum(value == 0 for value in node_state.values())
                infected_count = sum(value == 1 for value in node_state.values())
                recovered_count = sum(value == 2 for value in node_state.values())

                susceptible_count_t.append(susceptible_count)
                infected_count_t.append(infected_count)#1000步每步的感染人数数组
                recovered_count_t.append(recovered_count)

                f_infected_count_t = [value / N for value in infected_count_t]
                f_susceptible_count_t = [value / N for value in susceptible_count_t]
                f_recovered_count_t = [value / N for value in recovered_count_t]

            TT.append(f_infected_count_t)#一次迭代，1000步每步的感染比例数组'
            #print("f_infected_count_t:", f_infected_count_t)

        result_array = [sum(elements) / intera for elements in zip(*TT)]  # 将重复试验的数组元素取均值，得到新数组
        #print("result_array:", result_array)

        # 提取从第801个元素到最后一个元素之后的子数组，799表第800个元素
        sub_array = result_array[899:]
        # 计算子数组的平均值
        average = sum(sub_array) / len(sub_array)
        heatmap_data[delta_index, beta0_index] = average  # 取1000步时感染人数比例(-1表示最后一个时间步的数据)

# 保存 heatmap_data 到一个 NumPy 文件
np.save('heatmap_data.npy', heatmap_data)
plt.imshow(heatmap_data, cmap='jet', aspect='equal', origin='lower', interpolation='gaussian')
#plt.title('Smoothed Heatmap')
plt.xlabel('β(0)')
plt.ylabel('δ')
plt.colorbar()

# 自定义x轴刻度位置
x_tick_positions = [0, 0.2, 0.4, 0.6, 0.8, 1]
# 创建一个LinearLocator，以确保x轴刻度均匀分布
x_locator = LinearLocator(numticks=len(x_tick_positions))
plt.gca().xaxis.set_major_locator(x_locator)
# 设置x轴刻度标签
plt.gca().set_xticklabels([str(tick) for tick in x_tick_positions])

#自定义y轴刻度位置
y_tick_positions = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
# 创建一个LinearLocator，以确保y轴刻度均匀分布
y_locator = LinearLocator(numticks=len(y_tick_positions))
plt.gca().yaxis.set_major_locator(y_locator)
# 设置y轴刻度标签
plt.gca().set_yticklabels([str(tick) for tick in y_tick_positions])

plt.savefig('sir β&δ高斯热力图.pdf', format='pdf')
plt.pause(1)  # 在这里等待1秒钟
plt.show()


# 绘制原始热力图
#plt.subplot(1, 2, 2)
plt.imshow(heatmap_data, cmap='jet', aspect='auto',origin='lower')
plt.title('sir Original Heatmap')
# 添加颜色条
plt.colorbar()
# 添加坐标轴标签
plt.xlabel('β0')
#plt.ylabel('delta')
k3=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
k4=[0, 0.2, 0.4, 0.6, 0.8, 1]
# 设置x轴和y轴刻度 为beta0数组对应刻度
plt.yticks(np.arange(11), k3)
plt.xticks(np.arange(6), k4)
plt.savefig('sir β&δ原始热力图.pdf', format='pdf')
plt.show()

