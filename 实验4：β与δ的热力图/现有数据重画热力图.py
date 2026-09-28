import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import ndimage
from matplotlib.ticker import LinearLocator

heatmap_data = np.load('/论文1实验数据/实验4：β与δ的热力图/副本heatmap_data.npy')
plt.imshow(heatmap_data, cmap='jet', aspect='equal', origin='lower', interpolation='gaussian')
#plt.title('Smoothed Heatmap')
plt.xlabel('ζ')
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

plt.savefig('1.10重画副本β&δ高斯热力图.pdf', format='pdf')
plt.pause(1)  # 在这里等待1秒钟
plt.show()


# 绘制原始热力图
#plt.subplot(1, 2, 2)
plt.imshow(heatmap_data, cmap='jet', aspect='auto',origin='lower')
plt.title('Original Heatmap')
# 添加颜色条
plt.colorbar()
# 添加坐标轴标签
plt.xlabel('ζ')
#plt.ylabel('delta')

k3=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
k4=[0, 0.2, 0.4, 0.6, 0.8, 1]
# 设置x轴和y轴刻度 为beta0数组对应刻度
plt.yticks(np.arange(11), k3)
plt.xticks(np.arange(6), k4)
plt.savefig('1.10重画副本 β&δ原始热力图.pdf', format='pdf')
plt.show()