import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np


df = pd.read_table("data3.txt",skiprows=0)
fig, ax = plt.subplots()
density_profile = df.iloc[:,0]
print(len(density_profile))
#print(density_profile)
#x_range = np.arange(0,len(density_profile), step = 1)
x_range = np.arange(29,50, step = 1)
#print(x_range)
ax.plot(x_range, density_profile[29:50])
plt.grid()
plt.show()






