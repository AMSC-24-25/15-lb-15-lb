import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("lbm_results_3D.csv")

x, y, z, rho, ux, uy, uz = data['x'], data['y'], data['z'], data['rho'], data['ux'], data['uy'], data['uz']

# 选切片
z_middle = int(max(z) / 3)
mask = (z == z_middle)

# 提取切片数据
x_slice = x[mask]
y_slice = y[mask]
rho_slice = rho[mask]
speed_slice = np.sqrt(ux[mask]**2 + uy[mask]**2 + uz[mask]**2)

# 绘制速度大小切片
plt.figure(figsize=(8, 6))
plt.scatter(x_slice, y_slice, c=speed_slice, cmap='viridis')
plt.colorbar(label="Speed |u|")
plt.xlabel("X")
plt.ylabel("Y")
plt.savefig('Speed_3D.png')  
plt.close()

#运行 python print_image_3D.py 