from Spline import Spline
import matplotlib.pyplot as plt
import numpy as np
from time import time
import statistics

x = [1, 2, 3, 4, 5]
y = [1, 2, 3, 4, 5]

start_creation = time()
spline = Spline(x=x, y=y)
end_creation = time() - start_creation

ticks = np.linspace(0, 5, 50)
y_ticks = []

times_eval = []

for t in ticks:
    start_eval = time()
    v = spline.D0(t)
    end_eval = time() - start_eval
    times_eval.append(end_eval)
    y_ticks.append(v)

print('Tempo Creazione Spline:', end_creation)
print('Tempo Medio D0:', statistics.mean(times_eval))

# Tempo Creazione Spline: 0.0036389827728271484
# Tempo Medio D0: 5.545616149902344e-06


# plt.plot(x, y, 'o')
# plt.plot(ticks,y_ticks)
#
# plt.show()
