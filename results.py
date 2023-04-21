import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_excel('results.xlsx')

window = 1000
rms = df[['t1', 't2', 't3']].rolling(window).apply(lambda x: np.sqrt(np.mean(np.square(x))))

df_rms = rms.assign(Time=df.Time)

df.Time = df.Time / 60

fig, ax = plt.subplots()
ax.plot(df.Time, df.t1, label="1% Depth")
ax.plot(df.Time, df.t2, label="50% Depth")
ax.plot(df.Time, df.t3, label="90% Depth")

ax.set_xlabel('time [min]')
ax.set_ylabel('temperature [c]')
ax.set_title('Thermodynamic Siphon Experiment')
ax.legend()
plt.show()