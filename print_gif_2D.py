import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from PIL import Image  # 用于生成 GIF 动画

# 创建存储 PNG 文件和 GIF 动画的目录
output_directory = "./output_images/"
os.makedirs(output_directory, exist_ok=True)

# 设置迭代范围和步长
start_iter = 0
end_iter = 10000
step = 100

# 保存所有生成的 PNG 文件路径
png_files = []

# 固定网格尺寸
Nx = 100
Ny = 100

# 循环读取并绘制速度场
for iteration in range(start_iter, end_iter + 1, step):
    # 根据新的文件命名格式生成文件名
    filename = f"intermediate_step_{iteration}.csv"
    
    try:
        # 读取数据
        data = pd.read_csv(filename).dropna()  # 删除空行
        required_columns = ['rho', 'ux', 'uy']
        for column in required_columns:
            if column not in data.columns:
                raise ValueError(f"Column '{column}' not found in file '{filename}'.")

        # 检查数据大小是否匹配 100x100 网格
        total_points = len(data['rho'])
        if total_points != Nx * Ny:
            raise ValueError(f"Data size ({total_points}) does not match grid dimensions Nx={Nx}, Ny={Ny}.")

        # 提取并重塑数据
        rho = data['rho'].values.reshape(Ny, Nx)
        ux = data['ux'].values.reshape(Ny, Nx)
        uy = data['uy'].values.reshape(Ny, Nx)

        # 计算速度大小
        speed = np.sqrt(ux ** 2 + uy ** 2)

        # 绘制速度场
        plt.figure(figsize=(12, 8))  # 设置图像大小
        plt.imshow(speed, cmap='viridis', origin='lower')
        plt.colorbar(label='Speed', shrink=0.8)  # 调整色条大小
        plt.title(f'Speed Field (Iteration {iteration})')
        plt.axis('off')
        plt.tight_layout()  # 自动调整布局

        # 保存图片，文件名包含迭代次数
        png_path = os.path.join(output_directory, f"Speed_{iteration:04d}.png")
        plt.savefig(png_path, bbox_inches='tight', pad_inches=0.5, dpi=300)  # 增加边距和分辨率
        plt.close()

        # 将图片路径加入列表
        png_files.append(png_path)

        print(f"Speed field image for iteration {iteration} saved as '{png_path}'")

    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Skipping iteration {iteration}.")
    except ValueError as e:
        print(f"Error in file '{filename}': {e}. Skipping iteration {iteration}.")
    except Exception as e:
        print(f"Unexpected error in file '{filename}': {e}. Skipping iteration {iteration}.")

# 创建 GIF 动画
if png_files:
    gif_path = os.path.join(output_directory, "speed_evolution.gif")
    images = [Image.open(png) for png in png_files]
    images[0].save(
        gif_path,
        save_all=True,
        append_images=images[1:],
        duration=200,  # 每帧持续时间（毫秒）
        loop=0         # 无限循环
    )
    print(f"GIF animation saved as '{gif_path}'")
else:
    print("No PNG files were generated. GIF animation was not created.")
