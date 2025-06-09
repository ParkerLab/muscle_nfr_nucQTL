# args: 1. celltype 2. region 3.count matrix 4.output

# import pkgs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyreadr
import sys

# input
celltype = sys.argv[1]
region = sys.argv[2]
input = sys.argv[3]
output = sys.argv[4]
t1 = sys.argv[5]
t2 = sys.argv[6]
t3 = sys.argv[7]
t4 = sys.argv[8]
t5 = sys.argv[9]
t6 = sys.argv[10]
t7 = sys.argv[11]
t8 = sys.argv[12]

# readin
df_nucleo_x = pyreadr.read_r(input)[None]
# figure title
title = celltype + ' ' + region + ': Total peaks called=' + str(len(df_nucleo_x))

def feature_select(df, reads):
    xlist = range(df.shape[1])
    ylist = np.zeros(df.shape[1])
    m_color = np.zeros(df.shape)
    m_df = df.to_numpy()

    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            if m_df[i][j] >= reads:
                m_color[i][j] = 1

    color_summed = np.sum(m_color, axis=1)

    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            if (color_summed[i] >= xlist[j]):
                ylist[j] +=1
    return (xlist, ylist)

xlist, ylist1 = feature_select(df_nucleo_x, float(t1))
xlist, ylist2 = feature_select(df_nucleo_x, float(t2))
xlist, ylist3 = feature_select(df_nucleo_x, float(t3))
xlist, ylist4 = feature_select(df_nucleo_x, float(t4))
xlist, ylist5 = feature_select(df_nucleo_x, float(t5))
xlist, ylist10 = feature_select(df_nucleo_x, float(t6))
xlist, ylist15 = feature_select(df_nucleo_x, float(t7))
xlist, ylist20 = feature_select(df_nucleo_x, float(t8))

fig = plt.figure()
plt.plot(xlist, ylist1, label = t1)
plt.plot(xlist, ylist2, label = t2)
plt.plot(xlist, ylist3, label = t3)
plt.plot(xlist, ylist4, label = t4)
plt.plot(xlist, ylist5, label = t5)
plt.plot(xlist, ylist10, label = t6)
plt.plot(xlist, ylist15, label = t7)
plt.plot(xlist, ylist20, label = t8)
plt.xlabel('n')
plt.ylabel('features')
plt.legend(loc='right')
plt.title(title)
#plt.show()
# save fig
plt.tight_layout()
plt.savefig(output, dpi=300)