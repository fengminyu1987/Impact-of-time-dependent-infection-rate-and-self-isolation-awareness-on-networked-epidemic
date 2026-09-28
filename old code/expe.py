import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

# 定义WS网络参数
N = 1000  # 节点数量
K = 4  # 平均度
p = 0.5  # 重连概率

# 初始化WS网络
ws_graph = nx.watts_strogatz_graph(N, K, p)

# 初始化传染病模型，定义三个数组，元素个数分别为1, 0, N-1
infected_nodes = np.random.choice(N, size=1)  # 随机选择1个感染者，存储在下标为0的位置
recovered_nodes = []  # 初始状态下，没有康复者
susceptible_nodes = np.setdiff1d(np.arange(N), infected_nodes)  # 初始状态下，除一个感染者外都是易感者

beta0 = 0.1  # 初始传染率
mu = 0.5  # 康复率
alpha = 0.05  # 从r变为s的概率

# 模拟传播过程
total_time_steps = 10000
# 记录每个时间步的易感者、感染者和康复者数量
susceptible_count = []
infected_count = []
recovered_count = []
susceptible_count.append(len(susceptible_nodes))
infected_count.append(len(infected_nodes))
recovered_count.append(len(recovered_nodes))
for t in range(total_time_steps):  # 时间从0到n-1

    new_infected = []
    new_recovered = []
    new_susceptible = []

    for node in range(N):  # 0到n-1号节点
        if node in susceptible_nodes:

            # 与感染者相连的节点有一定概率成为新感染者
            neighbors = list(ws_graph.neighbors(node))
            infected_neighbor_count = sum(neighbor in infected_nodes for neighbor in neighbors)
            # for neighbor in neighbors即对于在neighbors中的元素neighbor，neighbor in infected_nodes即判断该元素
            # 是否也在感染节点列表，即是否为真（即布尔值），sum求和布尔值
            contact_probability = np.exp(-infected_neighbor_count / len(neighbors))

            #beta0 = 0.5  # 初始传染率
            sigma = 0.1  # 扩散率
            miu = 0  # 漂移率
            # 布朗运动模拟
            dt = 1  # 时间步长
            # 绘制传染率随时间变化的图像
            dB = np.sqrt(dt) * np.random.normal(size=total_time_steps)  # 取size=total_time_steps个数，Brownian motion increment,按正态概率正负无穷取值
            B = np.cumsum(dB)  # Brownian motion path，累加过后的数组
            beta = beta0 * np.exp((miu - 0.5 * sigma ** 2) * t + sigma * B)  # SIR模型中的传染率β的随机漂移

            gamma = contact_probability * beta[t]  # 计算概率 γ

            if np.random.random() < gamma:  # 生成一个随机数，如果随机数小于 gamma 的值，则将该节点添加到 new_infected 列表中。
                new_infected.append(node)

        elif node in infected_nodes:
            # 感染者以一定概率康复
            if np.random.random() < mu:
                new_recovered.append(node)

        elif node in recovered_nodes:
            # 康复者以一定概率变为易感者
            if np.random.random() < alpha:
                new_susceptible.append(node)

    # 更新感染者、康复者和易感者列表
    infected_nodes = np.concatenate((infected_nodes, new_infected))
    recovered_nodes = np.concatenate((recovered_nodes, new_recovered))
    susceptible_nodes = np.concatenate((susceptible_nodes, new_susceptible))

    # 移除旧的感染者、康复者和易感者节点
    infected_nodes = np.setdiff1d(infected_nodes, recovered_nodes)
    recovered_nodes = np.setdiff1d(recovered_nodes, susceptible_nodes)
    susceptible_nodes = np.setdiff1d(susceptible_nodes, infected_nodes)

    # 记录易感者、感染者和康复者数量
    susceptible_count.append(len(susceptible_nodes))
    infected_count.append(len(infected_nodes))
    recovered_count.append(len(recovered_nodes))

# 绘制易感者、感染者和康复者数量随时间的变化曲线
f_infected_count = [value/N for value in infected_count]
plt.plot(range(total_time_steps+1), f_infected_count,color='red', label='Fi with sid and di')

#plt.plot(range(total_time_steps+1), susceptible_count, label='susceptible')
#plt.plot(range(total_time_steps+1), infected_count, label='Infected')
#plt.plot(range(total_time_steps+1), recovered_count, label='Recovered')

plt.xlabel('Time Steps')
plt.ylabel('Number of Nodes')
plt.ylim(0,1)
plt.legend()# 添加图例，即用来解释曲线代表的东西
#plt.xscale('log')
plt.title('expe, Experimental group vs Control group')
plt.grid()  # 添加网格线
plt.show()

f_infected_count4 = f_infected_count










