from Spline import Spline
# import matplotlib.pyplot as plt
import numpy as np
from time import time
import statistics

x = [1, 2, 3, 4, 5]
y = [1, 2, 3, 4, 5]

start_creation = time()
spline = Spline(x=x, y=y, verbose=False)
end_creation = time() - start_creation

ticks = np.linspace(0, 5, 5000)
y_ticks = []
y1_ticks = []
y2_ticks = []

times_eval = []

for t in ticks:
    start_eval = time()
    v = spline.evaluate(t)
    end_eval = time() - start_eval
    times_eval.append(end_eval)
    y_ticks.append(v)
    y1_ticks.append(spline.evaluate(t, der=1))
    y2_ticks.append(spline.evaluate(t, der=2))

print('Tempo Creazione Spline:', end_creation)
print('Tempo Medio D0:', statistics.mean(times_eval))

# Senza flag compilatore
# Tempo Creazione Spline: 0.0036389827728271484
# Tempo Medio D0: 5.545616149902344e-06

# Con flag
# Tempo Creazione Spline: 0.0013260841369628906
# Tempo Medio D0: 4.644393920898438e-06


# plt.plot(x, y, 'o')
# plt.plot(ticks, y_ticks)
# plt.plot(ticks, y1_ticks)
# plt.plot(ticks, y2_ticks)
#
# plt.show()
