from SplinePoliMi import Spline
import matplotlib.pyplot as plt
import numpy as np
from time import time


x = [6, 2, 3, 4, 5, 5, 6]
y = [1, 4, 16, 4, 1, 2, 1]

start_creation = time()
spline = Spline(x=x, y=y, verbose=False)
end_creation = time() - start_creation

ticks = np.linspace(0, 5, 5000)


start_eval = time()
y_ticks = spline.evaluate(ticks, der=0)
y1_ticks = spline.evaluate(ticks, der=1)
y2_ticks = spline.evaluate(ticks, der=2)
end_eval = time() - start_eval


print('Tempo Creazione Spline:', end_creation)
print('Tempo Medio D0:', end_eval)

# Senza flag compilatore
# Tempo Creazione Spline: 0.0036389827728271484
# Tempo Medio D0: 5.545616149902344e-06

# Con flag
# Tempo Creazione Spline: 0.0013260841369628906
# Tempo Medio D0: 4.644393920898438e-06


plt.plot(x, y, 'o')
plt.plot(ticks, y_ticks)
plt.plot(ticks, y1_ticks)
plt.plot(ticks, y2_ticks)

plt.show()
