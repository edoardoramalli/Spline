from SplinePoliMi import Spline
import matplotlib.pyplot as plt
import numpy as np
from time import time


# x = [6, 2, 3, 4, 5, 5, 6]
# y = [1, 4, 16, 4, 1, 2, 1]


x = [560.0, 590.0, 620.0, 650.0, 680.0, 710.0, 740.0, 770.0, 800.0, 870.0, 945.0, 1030.0]
y = [0.000486207, 0.000262069, 0.000144828, 0.000131034, 0.000168966, 0.000224138, 0.000272414, 0.00032069, 0.000217241, 2.75862e-05, 3.4483e-06, 3.4483e-06]


start_creation = time()
spline = Spline(x=x, y=y)

end_creation = time() - start_creation

ticks = np.linspace(min(x), max(x), 250)


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
