import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skimage.measure import marching_cubes

data = pd.read_csv("result/lbm_results_3D.csv")

Nx, Ny, Nz = 100, 100, 100
rho = data['rho'].values.reshape(Nz, Ny, Nx)
ux = data['ux'].values.reshape(Nz, Ny, Nx)
uy = data['uy'].values.reshape(Nz, Ny, Nx)
uz = data['uz'].values.reshape(Nz, Ny, Nx)
speed = np.sqrt(ux**2 + uy**2 + uz**2)


iso_value = np.median(speed)

verts, faces, normals, values = marching_cubes(speed, level=iso_value)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')

mesh = ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2],
                       cmap='viridis', lw=1)
fig.colorbar(mesh, label='Speed')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Speed Isosurface')
os.makedirs("picture", exist_ok=True)
plt.savefig("picture/3D_Speed_Isosurface.png")
plt.close()
