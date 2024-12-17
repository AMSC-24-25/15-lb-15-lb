import pandas as pd

# 读取原始 CSV 文件
data = pd.read_csv("lbm_results.csv")

# 获取数据总条数
total_rows = len(data)

# 分割数据并保存为多个文件
for end in range(100, total_rows + 1, 100):
    # 提取从第一条到第 `end` 条数据
    subset = data.iloc[:end]
    
    # 保存到新的 CSV 文件，文件名包含当前迭代步数
    subset.to_csv(f"lbm_results_{end}.csv", index=False)
    print(f"File 'lbm_results_{end}.csv' saved.")
