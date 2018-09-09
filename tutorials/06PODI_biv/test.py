import numpy as np
import matplotlib.pyplot as plt


# mu_all = np.array([0.3333333333, 0.2, 0.1428571429, 0.1111111111, 0.0909090909, 0.0769230769, 0.0666666667,
#                0.0588235294, 0.0526315789, 0.0476190476, 0.0434782609, 0.04, 0.037037037, 0.0344827586, 0.0322580645,
#                0.0303030303, 0.0285714286, 0.027027027, 0.0256410256, 0.0243902439, 0.023255814, 0.0222222222, 0.0212765957,
#                0.0204081633, 0.0196078431, 0.0188679245, 0.0181818182, 0.0175438596, 0.0169491525, 0.0163934426, 0.0158730159,
#                0.0153846154, 0.0149253731, 0.0144927536, 0.014084507, 0.0136986301, 0.0133333333, 0.012987013, 0.0126582278,
#                0.012345679, 0.0120481928, 0.0117647059, 0.0114942529, 0.0112359551, 0.010989011, 0.0107526882, 0.0105263158,
#                0.0103092784, 0.0101010101])
mu_all = np.array([0.0322580645, 0.0303030303, 0.0285714286, 0.027027027, 0.0256410256, 0.0243902439, 0.023255814, 0.0222222222, 0.0212765957,
               0.0204081633, 0.0196078431, 0.0188679245, 0.0181818182, 0.0175438596, 0.0169491525, 0.0163934426, 0.0158730159,
               0.0153846154, 0.0149253731, 0.0144927536, 0.014084507, 0.0136986301, 0.0133333333, 0.012987013, 0.0126582278,
               0.012345679, 0.0120481928, 0.0117647059, 0.0114942529, 0.0112359551, 0.010989011, 0.0107526882, 0.0105263158,
               0.0103092784, 0.0101010101])
time_all = np.arange(10,20.01,0.5)

print(mu_all.size)

mu_given = 5./mu_all[0::2]
mu_interp = 5./mu_all[1::2]
print(mu_given.size)
print(mu_interp.size)

time_given = time_all[0::2]
time_interp = time_all[1::2]

# fig, ax = plt.subplots()
# fig = plt.figure(figsize=(80, 60), dpi=100, facecolor='w', edgecolor='k')
fig = plt.figure(figsize=(100, 80), dpi=100, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)

time_given, mu_given = np.meshgrid(time_given, mu_given)
time_interp, mu_interp = np.meshgrid(time_interp, mu_interp)
ax.scatter(time_given, mu_given, c='red', label='Given snapshots', alpha=0.3, edgecolors='none')
ax.scatter(time_interp, mu_interp, c='blue', label='Required snapshots', alpha=0.3, edgecolors='none')

title = 'Illustrative grid for the given and required snapshots'
plt.title(title)
plt.xlabel('TIME')
plt.ylabel('REYNOLDS NUMBER')

# ax.legend(loc=2)
ax.legend(loc=9, bbox_to_anchor=(0.5, -0.075), ncol=2)

# Major ticks every 20, minor ticks every 5
major_ticks = np.arange(10, 20.1, 1)
minor_ticks = np.arange(10.5, 20, 1)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
major_ticks = 5./mu_all[0::2]
minor_ticks = 5./mu_all[1::2]
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)

# And a corresponding grid
ax.grid(which='both')

# Or if you want different settings for the grids:
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)

ax.grid(True, linestyle='dotted')

plt.show()