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

fig, axs = plt.subplots(3, figsize=(10, 15))

versions = []
modes = []
for h in y:
    versions.append(h['name'])
    modes += h['modes'].keys()
modes = sorted(set(modes))
# print(modes)

plt.rcParams['image.cmap']='jet'

axs[0].grid(c='0.9')
axs[0].set_axisbelow(True)
axs[0].set_title('Run time, s (lower better)', fontstyle='italic')

axs[1].grid(c='0.9')
axs[1].set_axisbelow(True)
axs[1].set_title('Incremental compile time, s (lower better)', fontstyle='italic')

axs[2].grid(c='0.9')
axs[2].set_axisbelow(True)
axs[2].set_title('Coefficient: run time divided by release time, 0..1 (higher better)', fontstyle='italic')

for mode in modes:
    keys = []
    times = []
    comp_time = []
    comp_time2 = []
    cf = []
    for h in y:
        version = h['name']
        mode_h = h['modes']
        if mode in mode_h:
            keys.append(version)
            times.append(mode_h[mode]['time'])
            comp_time.append(mode_h[mode]['compile1'])
            comp_time2.append(mode_h[mode]['compile2'])
            cf.append(mode_h[mode]['cf'])
    axs[0].plot(keys, times)
    axs[0].scatter(keys, times, label=mode, s=100, marker='D')

    axs[1].plot(keys, comp_time2)
    axs[1].scatter(keys, comp_time2, label=mode, s=100, marker='D')

    if "release" not in mode:
        axs[2].plot(keys, cf)
        axs[2].scatter(keys, cf, label=mode, s=100, marker='D')

plt.grid(color='0.95')
axs[0].legend(loc='center left')
axs[1].legend(loc='center left')
axs[2].legend(loc='center left')
plt.savefig('releases.png')
