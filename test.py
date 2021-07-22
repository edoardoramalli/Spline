from Spline import Spline
import matplotlib.pyplot as plt
import numpy as np

x = [1, 2, 3, 4, 5]
y = [1, 2, 3, 4, 5]

spline = Spline(x=x, y=y)

ticks = np.linspace(0, 5, 50)

plt.plot(x, y, 'o')
plt.plot(ticks, [spline.D0(t) for t in ticks])

plt.show()
