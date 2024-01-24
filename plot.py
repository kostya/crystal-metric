# pip3 install matplotlib
import matplotlib.pyplot as plt
# pip3 install pyyaml
import yaml

with open("history.yml", "r") as stream:
    try:
        y = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

fig, ax = plt.subplots()

for key in y.keys():
    l = y[key]
    releases = []
    values = []
    for k, v in l:
        releases.append(k)
        values.append(v)
    ax.plot(releases, values, label=key)
    ax.scatter(releases, values)
    
plt.legend()
#plt.show()
plt.savefig('releases.png')
