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

intera = 50
delta_array = [0, 1, 2, 3, 4, 5, 6,7,8]
final_result_array = []
for delta in delta_array:
    print(delta)
    print(str(delta/100*100)+"%")
    TT = [] #针对每个beta值，收集 五次迭代 后的数组
    for _ in range(intera):
        node_state = {}
        # 0: S; 1: I; 2: R.
        for node in ws_graph.nodes:
            node_state[node] = 0
        initial_infected_nodes = random.sample(list(ws_graph.nodes), 10)  # 随机选择100个节点
        for node in initial_infected_nodes:
            node_state[node] = 1

        beta0 = 0.5  # S->I
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
                           xy_ratio = 1- (np.exp(- delta * infected_neis_count/neigh_num)*beta0)
                           total_xy_ratio *= xy_ratio #将每个感染邻居的（1-cβ）值相乘
                        else:
                            continue
                    probai = 1-total_xy_ratio

                    if 0 < random.uniform(0, 1) < probai:
                        node_state[node]=1

                elif node_state[node]==1:
                    if 0<random.uniform(0,1)<gamma:
                        node_state[node]= 2

                elif node_state[node] == 2:
                    node_state[node] = 2

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
    plt.plot(range(total_time_steps+1), final_result_array[i])


plt.xlabel('Time')
plt.ylabel('Infection Density')
#plt.xscale('log')
plt.xlim(0, 101)
plt.ylim(0.01, 1)  # 设置 y 轴的显示范围，从 0.01 到 0.6

#plt.yticks(np.arange(0, 0.61, 0.1))
#plt.legend() # 添加图例，即用来解释曲线代表的东西
#plt.title(f'Time Evolution')
plt.grid()  # 添加网格线

plt.savefig('eg实验52——无dy SIR不同δ下实验.pdf', format='pdf')
plt.show()


f_infected_count_exp = f_infected_count_t

import csv#保存每个δ对应感染密度最大值

# 文件路径
file_path = 'eg-sir_无Dir_max_infection_values.csv'

# 将最大值写入CSV文件
with open(file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Delta', 'Max Infection'])
    for i, result_array in enumerate(final_result_array):
        max_infection = max(result_array)
        writer.writerow([delta_array[i], max_infection])

print("最大感染曲线值已保存到:", file_path)













