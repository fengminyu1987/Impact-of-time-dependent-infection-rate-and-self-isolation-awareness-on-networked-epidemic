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
beta0_array = [0.1, 0.3, 0.5, 0.7, 0.9]
final_result_array = []
for beta0 in beta0_array:
    print(beta0)
    TT = []
    for _ in range(intera):
        node_state = {}
        # 0: S; 1: I; 2: R.
        for node in ws_graph.nodes:
            node_state[node] = 0
        initial_infected_nodes = random.sample(list(ws_graph.nodes), 10)  # 随机选择1%个节点
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

        for t in range(total_time_steps):
            for node in ws_graph.nodes:
                if node_state[node]==0:
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
                           xy_ratio = 1- (np.exp(-infected_neis_count/neigh_num)*beta0)
                           total_xy_ratio *= xy_ratio #将每个感染邻居的（1-cβ）值相乘
                        else:
                            continue
                    probai = 1-total_xy_ratio

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
            f_infected_count_t = [value / N for value in infected_count_t]
            f_susceptible_count_t = [value / N for value in susceptible_count_t]
            f_recovered_count_t = [value / N for value in recovered_count_t]

        TT.append(f_infected_count_t)

    result_array = [sum(elements)/intera for elements in zip(*TT)]
    final_result_array.append(result_array) #元素个数等于beta数组个数
for i in range(len(final_result_array)):
    plt.plot(range(total_time_steps+1), final_result_array[i], label="β={}".format(beta0_array[i]))

plt.xlabel('Time')
plt.ylabel('Infection Density')
plt.xscale('log')
plt.xlim(0, 1001)
plt.yticks(np.arange(0, 0.81, 0.1))
plt.legend() # 添加图例，即用来解释曲线代表的东西
plt.title('Time Evolution')
plt.grid()  # 添加网格线

plt.savefig('感染率上升_无动感有自保.pdf', format='pdf')
plt.show()



f_infected_count_exp = f_infected_count_t








