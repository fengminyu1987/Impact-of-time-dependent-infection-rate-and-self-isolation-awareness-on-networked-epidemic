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

intera = 10
beta0_array = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
final_result_array = []
for beta0 in beta0_array:
    print("beta0:", beta0)
    TT = [] #针对每个beta值，收集 五次迭代 后的数组
    for _ in range(intera):
        node_state = {}
        # 0: S; 1: I; 2: R.
        for node in ws_graph.nodes:
            node_state[node] = 0
        #node_state[random.choice(list(ws_graph.nodes))] = 1
        initial_infected_nodes = random.sample(list(ws_graph.nodes), 10)  # 随机选择10个节点
        for node in initial_infected_nodes:
            node_state[node] = 1

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

        # 动态感染率
        sigma = 0.02  # 扩散率
        miu = -0.01  # 漂移率
        # 布朗运动模拟
        dt = 1  # 时间步长
        # 绘制传染率随时间变化的图像
        dB = np.sqrt(dt) * np.random.normal(size=total_time_steps)  # 取size=total_time_steps个数，Brownian motion increment,按正态概率正负无穷取值
        B = np.cumsum(dB)  # Brownian motion path，累加过后的数组
        beta = [0]*total_time_steps
        for t in range(total_time_steps):
            beta[t] = beta0 * np.exp((miu - 0.5 * sigma ** 2) * t + sigma * B[t])   # SIR模型中的传染率β的随机漂移

        for t in range(total_time_steps):
            for node in ws_graph.nodes:
                if node_state[node]==0:
                    #自我隔离强度
                    neis=list(ws_graph.neighbors(node))
                    total_xy_ratio=1
                    for neighbor1 in neis:  #遍历一代邻居
                        if node_state[neighbor1]==1:
                            # 自我隔离强度
                            neis = list(ws_graph.neighbors(node))
                            total_xy_ratio = 1
                            for neighbor1 in neis:  # 遍历一代邻居
                                if node_state[neighbor1] == 1:  # 若一代邻居不正常，计数 2代正常和感染人数，1-（1-beta）^n
                                    neis2 = list(ws_graph.neighbors(neighbor1))
                                    neigh_num = len(neis2)  # 2代邻居数
                                    infected_neis_count = 0
                                    for neighbor2 in neis2:  # 遍历2代邻居，找不正常数
                                        if node_state[neighbor2] == 1:
                                            infected_neis_count += 1  # 2代感染数
                                        else:
                                            continue
                                    xy_ratio = 1 - (np.exp(-infected_neis_count / neigh_num) * beta[t])  # 没被2代感染概率
                                    total_xy_ratio *= xy_ratio  # 总的没被感染率，每（1-cβ）值相乘
                                    print("打印 beta[t]")
                                    print(beta[t])
                                else:
                                    continue  # 若某一代邻居正常，则不用管，判断下个一代邻居
                    probai = 1 - total_xy_ratio  # 若一代邻居都正常，probai=1-1=0

                    if 0 < random.uniform(0, 1) < probai:
                        node_state[node]=1

                elif node_state[node]==1:
                    if 0<random.uniform(0,1)<gamma:
                        node_state[node]=2

                elif node_state[node]==2:
                    if 0<random.uniform(0,1)<alpha:
                        node_state[node]=0

            #print(f"计算{t}时刻sir态个体数目并保存在数组中")
            susceptible_count = sum(value == 0 for value in node_state.values())
            infected_count = sum(value == 1 for value in node_state.values())
            recovered_count = sum(value == 2 for value in node_state.values())

            susceptible_count_t.append(susceptible_count)
            infected_count_t.append(infected_count)
            recovered_count_t.append(recovered_count)

        # 绘制易感者、感染者和康复者数量随时间的变化曲线
            f_infected_count_t = [value / N for value in infected_count_t] #累加元素的数组
            f_susceptible_count_t = [value / N for value in susceptible_count_t]
            f_recovered_count_t = [value / N for value in recovered_count_t]

        TT.append(f_infected_count_t)
    result_array = [sum(elements)/intera for elements in zip(*TT)]#将重复试验得到的数组的元素取均值，得到新数组
    final_result_array.append(result_array) #元素个数等于beta数组个数
for i in range(len(final_result_array)):
    plt.plot(range(total_time_steps+1), final_result_array[i], label="β={}".format(beta0_array[i]))

plt.xlabel('Time')
plt.ylabel('Infection Density')
plt.xscale('log')
plt.yticks(np.arange(0, 0.81, 0.1))
plt.xlim(0, 1001)
plt.legend() # 添加图例，即用来解释曲线代表的东西
plt.title('Time Evolution')
plt.grid()  # 添加网格线

plt.savefig('11月动感加自保.pdf', format='pdf')
plt.show()

f_infected_count_exp = f_infected_count_t








