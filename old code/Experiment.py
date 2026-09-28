import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
# 定义WS网络参数
N = 1000  # 节点数量
K = 4  # 平均度
p = 0.5  # 重连概率
# 初始化WS网络
ws_graph = nx.watts_strogatz_graph(N, K, p)

# 模拟传播过程
total_time_steps = 1000

node_state={}
#0: S; 1: I; 2: R.
for node in ws_graph.nodes:
    node_state[node]=0
node_state[random.choice(list(ws_graph.nodes))]=1

beta0=0.5#S->I
gamma=0.4#I->R
alpha=0.5#R->S

susceptible_count_t = []
infected_count_t = []
recovered_count_t = []
susceptible_count = sum(value == 0 for value in node_state.values())
infected_count = sum(value == 1 for value in node_state.values())
recovered_count = sum(value == 2 for value in node_state.values())
#计算初始态sir个数
susceptible_count_t.append(susceptible_count)
infected_count_t.append(infected_count)
recovered_count_t.append(recovered_count)

for t in range(total_time_steps):
    for node in ws_graph.nodes:
        if node_state[node]==0:
            #动态感染率
            sigma = 0.001  # 扩散率
            miu = 0  # 漂移率
            # 布朗运动模拟
            dt = 1  # 时间步长
            # 绘制传染率随时间变化的图像
            dB = np.sqrt(dt) * np.random.normal(size=total_time_steps)  # 取size=total_time_steps个数，Brownian motion increment,按正态概率正负无穷取值
            B = np.cumsum(dB)  # Brownian motion path，累加过后的数组
            beta = beta0 * np.exp((miu - 0.5 * sigma ** 2) * t + sigma * B)  # SIR模型中的传染率β的随机漂移

            #自我隔离强度
            neis=list(ws_graph.neighbors(node))
            total_xy_ratio=1
            for neighbor1 in neis:  #遍历一代邻居
                if node_state[neighbor1]==1:
                   neigh_num = 0
                   infected_neis_count = 0
                   for neighbor2 in ws_graph.neighbors(neighbor1): #遍历一代邻居的邻居
                        neigh_num+=1
                        if node_state[neighbor2]==1:
                           infected_neis_count += 1
                   xy_ratio = 1- (np.exp(-infected_neis_count/neigh_num)*beta[t])
                   total_xy_ratio *= xy_ratio #将每个感染邻居的（1-cβ）值相乘
                else:
                    continue

            probai = 1-total_xy_ratio

            if 0 < random.uniform(0, 1) < probai:
                node_state[node]=1

        if node_state[node]==1:
            if 0<random.uniform(0,1)<gamma:
                node_state[node]=2

        if node_state[node]==2:
            if 0<random.uniform(0,1)<alpha:
                node_state[node]=0

    print(f"计算{t}时刻sir态个体数目并保存在数组中")
    susceptible_count = sum(value == 0 for value in node_state.values())
    infected_count = sum(value == 1 for value in node_state.values())
    recovered_count = sum(value == 2 for value in node_state.values())
    susceptible_count_t.append(susceptible_count)
    infected_count_t.append(infected_count)
    recovered_count_t.append(recovered_count)

# 绘制易感者、感染者和康复者数量随时间的变化曲线
f_infected_count_t = [value / N for value in infected_count_t]
f_susceptible_count_t = [value / N for value in susceptible_count_t]
f_recovered_count_t = [value / N for value in recovered_count_t]

plt.plot(range(total_time_steps + 1), f_infected_count_t, color='red', label='Fi with sid and di')

'''
plt.plot(range(total_time_steps + 1), f_susceptible_count_t, color='green', label='Fs with sid and di')
plt.plot(range(total_time_steps + 1), f_recovered_count_t, color='blue', label='Fr with sid and di')

plt.plot(range(total_time_steps+1), susceptible_count_t, label='susceptible')
plt.plot(range(total_time_steps+1), infected_count_t, label='Infected')
plt.plot(range(total_time_steps+1), recovered_count_t, label='Recovered')
'''
plt.xlabel('Time Steps')
plt.ylabel('Number of Nodes')
plt.legend() # 添加图例，即用来解释曲线代表的东西
#plt.xscale('log')
plt.title(f'Experimental group with beta0={beta0},gamma={gamma},alpha={alpha},sigam{sigma}')
plt.grid()  # 添加网格线
plt.show()

f_infected_count_exp = f_infected_count_t










