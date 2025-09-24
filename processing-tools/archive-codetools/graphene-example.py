import numpy as np
import matplotlib.pyplot as plt

a1 = np.array([1,0])
a2 = np.array([1/2, np.sqrt(3)/2])



p = a1+a2
p2 = -a1-a2

x = [1, a2[0], p[0], p2[0]]
y = [0, a2[1], p[1], p2[1]]

plt.figure(figsize=(6,6), dpi=200)
plt.scatter(x, y)
plt.scatter([0, 2/3],[0, 1/3], label='atoms')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.legend()
plt.show()