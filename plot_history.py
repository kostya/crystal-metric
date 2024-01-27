# pip3 install matplotlib
import matplotlib.pyplot as plt
# pip3 install pyyaml
import yaml
from packaging import version

with open("history.yml", "r") as stream:
    try:
        y = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

fig, axs = plt.subplots(2, figsize=(10, 10))

versions = sorted(y.keys(), key=lambda v: version.parse(v))
modes = []
for version in y.keys():
    modes += y[version]
modes = sorted(set(modes))
# print(modes)

plt.rcParams['image.cmap']='jet'

axs[0].grid(c='0.9')
axs[0].set_axisbelow(True)
axs[0].set_title('Run time, s', fontstyle='italic')

axs[1].grid(c='0.9')
axs[1].set_axisbelow(True)
axs[1].set_title('Incremental compile time, s', fontstyle='italic')

for mode in modes:
    keys = []
    times = []
    comp_time = []
    comp_time2 = []
    for version in versions:
        if mode in y[version]:
            keys.append(version)
            times.append(y[version][mode]['time'])
            comp_time.append(y[version][mode]['compile1'])
            comp_time2.append(y[version][mode]['compile2'])
    axs[0].plot(keys, times)
    axs[0].scatter(keys, times, label=mode, s=100, marker='D')
    axs[1].plot(keys, comp_time2)
    axs[1].scatter(keys, comp_time2, label=mode, s=100, marker='D')

plt.grid(color='0.95')
axs[0].legend(loc='center left')
axs[1].legend(loc='center left')
plt.savefig('releases.png')
