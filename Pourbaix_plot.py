# coding=utf-8
#  Copyright (c) 2025. By ZYF
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy import constants
from mpl_toolkits.axes_grid1 import make_axes_locatable


#Data
kb   = constants.Boltzmann
T    = 298.15
e    = constants.e
ln10 = np.log(10)
epsilon = 2
C_H = 25
epsilon_0 = constants.epsilon_0
U_PZC = 0.3

G_H2O= -14.22
G_H2 = -6.8
mu = 0.8
alpha = 0.8

def G_sp(G_bare, G_tot, pH, U_RHE, m, n, classic=True):
    G_sp = G_bare + m * G_H2O - G_tot - (2 * m - n) * (G_H2 / 2 - U_RHE)
    if not classic:
        E = (U_RHE - kb * T * ln10 * pH - U_PZC) / (epsilon * epsilon_0)
        G_sp = G_sp + mu * E + alpha / 2 * E**2
    return G_sp

def select(x, y):
    G_bare = -279.53  #bare site energy
    G_tot = np.array([-282.11, G_bare, -294.44, -289.82, -285.48])  #total energy
    m = np.array([0, 0, 2, 1, 1])  #Count of O atom
    n = np.array([1, 0, 1, 1, 0])  #Count of H atom
    z = np.zeros_like(G_tot)
    for i in range(z.shape[0]):
        z[i] = G_sp(G_bare, G_tot[i], x, y, m[i], n[i], classic=True)
    #print(x, y, z)
    return np.argmin(-z)

def create_orr_pourbaix_diagram():
    # 创建网格
    DPI = 300
    pH = np.linspace(-1, 14, DPI)
    U = np.linspace(-2, 2, DPI)  # ORR相关电势范围
    X, Y = np.meshgrid(pH, U)

    # 创建数据矩阵
    Z = np.zeros_like(X)
    for x in range(DPI):
        for y in range(DPI):
            Z[y, x] = select(pH[x], U[y])


    # 创建颜色映射
    # 创建自定义颜色映射
    colors = [
        "#E78AC3",  # 樱花粉
        "#8DA0CB",  # 淡紫色
        "#66C2A5",  # 海绿色
        "#FC8D62",  # 珊瑚橙
        "#FFD92F"  # 柔黄色
    ]
    custom_cmap = ListedColormap(colors)

    # 创建离散化颜色映射
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
    norm = BoundaryNorm(bounds, custom_cmap.N)

    # 创建图形
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(right=0.85)

    # 绘制Pourbaix图
    pourbaix = ax.pcolormesh(X, Y, Z,
                             cmap=custom_cmap,
                             norm=norm,
                             shading='auto',
                             alpha=0.85)

    # 添加水稳定区边界
    ax.plot(pH, -0.059 * pH, 'k-', linewidth=2.5, label='HER (0V)')
    ax.plot(pH, 1.23 - 0.059 * pH, 'k--', linewidth=2.5, label='OER (1.23V)')

    # 设置坐标轴范围和标签
    ax.set_xlim(0, 7)
    ax.set_ylim(-2, 2)
    ax.set_xlabel('pH', fontsize=14, fontweight='bold')
    ax.set_ylabel('Potential (V vs. RHE)', fontsize=14, fontweight='bold')
    ax.set_title('Pourbaix Diagram with Color Legend', fontsize=16, fontweight='bold')

    # 添加网格
    ax.grid(True, linestyle='--', alpha=0.3)

    # ========== 添加右侧色谱图例 ==========
    # 创建图例坐标轴
    # 创建可分割的坐标轴
    divider = make_axes_locatable(ax)

    # 在右侧添加新的坐标轴（5%宽度，0.3英寸间距）
    cax = divider.append_axes("right", size="5%", pad=0.3)

    # 创建颜色条
    cbar = fig.colorbar(pourbaix, cax=cax)

    # 设置图例标签
    region_labels = [
        '*H',
        'slab',
        '*OOH',
        '*OH',
        '*O'
    ]

    # 设置图例刻度位置（在每个颜色块中间）
    cbar.set_ticks([0, 1, 2, 3, 4])
    cbar.ax.set_yticklabels(region_labels, fontsize=10)
    cbar.set_label('Stable Phases', fontsize=12, fontweight='bold')

    # 美化图例
    cax.tick_params(size=0)  # 隐藏刻度线
    cbar.outline.set_visible(False)  # 隐藏边框

    # 添加关键标注
    #ax.text(2, -1.2, 'H₂ Evolution', fontsize=10, rotation=-12, color='blue')
    #ax.text(2, 1.5, 'O₂ Evolution', fontsize=10, rotation=-12, color='red')
    #ax.text(10, 0.2, 'Water Stability Region', fontsize=12, bbox=dict(facecolor='white', alpha=0.8))

    # 添加区域标注
    #ax.text(2, -1.5, 'Metal Stability', fontsize=10, color='white', weight='bold')
    #ax.text(6, 0.0, 'Ion Stability', fontsize=10, color='black', weight='bold')
    #ax.text(10, 0.0, 'Hydroxide', fontsize=10, color='black', weight='bold')
    #ax.text(7, 1.2, 'Oxide Formation', fontsize=10, color='black', weight='bold')

    plt.tight_layout()
    plt.savefig('orr_pourbaix.png', dpi=300)
    plt.show()


# 生成ORR Pourbaix图
create_orr_pourbaix_diagram()