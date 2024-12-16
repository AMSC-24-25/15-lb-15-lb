import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("lbm_results_3D.csv")

Nx, Ny, Nz = 100, 100, 100

# 提取数据
rho = data['rho'].values.reshape(Nz, Ny, Nx)   # Density
ux = data['ux'].values.reshape(Nz, Ny, Nx)     # x 
uy = data['uy'].values.reshape(Nz, Ny, Nx)     # y 
uz = data['uz'].values.reshape(Nz, Ny, Nx)     # z 

# 计算速度大小 
speed = np.sqrt(ux**2 + uy**2 + uz**2)

# 选择切片 
z_slice = 50
rho_slice = rho[z_slice, :, :]     
speed_slice = speed[z_slice, :, :] 

# Density
plt.figure(figsize=(8, 6))
plt.imshow(rho_slice, cmap='viridis', origin='lower')
plt.colorbar(label='Density')
tile_Density='Density Field (3D) '+'z_slice = '+str(z_slice)
plt.title(tile_Density)
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig('Density_Field.png')
plt.close()

# Speed
plt.figure(figsize=(8, 6))
plt.imshow(speed_slice, cmap='plasma', origin='lower')
plt.colorbar(label='Speed')
tile_speed='Speed Field (3D) '+'z_slice = '+str(z_slice)
plt.title(tile_speed)
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig('Speed_Field.png')
plt.close()

