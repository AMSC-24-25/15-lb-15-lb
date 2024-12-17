import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("lbm_results.csv")

Nx, Ny = 100, 100
rho = data['rho'].values.reshape(Ny, Nx)
ux = data['ux'].values.reshape(Ny, Nx)
uy = data['uy'].values.reshape(Ny, Nx)
speed = np.sqrt(ux ** 2 + uy ** 2)

plt.imshow(speed, cmap='viridis', origin='lower')
plt.colorbar(label='Speed')
plt.title('Speed Field')
plt.savefig('Speed.png')  
plt.close()
