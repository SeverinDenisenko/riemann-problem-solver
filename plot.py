from matplotlib import pyplot as plt
import pandas as pd

df1 = pd.read_csv("first.dat", sep="\s+")
df2 = pd.read_csv("second.dat", sep="\s+")
df3 = pd.read_csv("third.dat", sep="\s+")
df4 = pd.read_csv("sod.dat", sep="\s+")

title1 = "{ro, v, p}L = {1, 0, 3}, {ro, v, p}R = {1, 0, 1}"
title2 = "{ro, v, p}L = {1, 1, 3}, {ro, v, p}R = {1,−1, 1}"
title3 = "{ro, v, p}L = {1,−0.1, 1}, {ro, v, p}R = {1, 0.2, 1}"
title4 = "{ro, v, p}L = {1, 0, 1}, {ro, v, p}R = {0.125, 0, 0.1}"

for (df, title) in ((df1, title1), (df2, title2), (df3, title3), (df4, title4)):
    fig, (v, ro, p) = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    fig.suptitle(title)
    v.set_title("velocity")
    v.plot(df["x"], df["v"])
    ro.set_title("density")
    ro.plot(df["x"], df["ro"])
    p.set_title("pressure")
    p.plot(df["x"], df["p"])
    plt.tight_layout()
    plt.show()
