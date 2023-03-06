import numpy as np 
import matplotlib.pyplot as plt

plt.style.use('bmh')

title = "Performance plot"
flags = ""
unity = "flops/cycle"

if unity == "" or title == "":
       print("Missing something")
       exit()


# Plot values
t = np.array([2**i for i in range(4,24)])

values = np.array([133.000000,
       73.000000,
       261.000000,
       137.000000,
       517.000000,
       265.000000,
       1030.000000,
       521.000000,
       2061.500000,
       1034.000000,
       4114.000000,
       2067.000000,
       8242.500000,
       4143.000000,
       16560.500000,
       8355.500000,
       32878.000000,
       16437.000000,
       65637.000000,
       32830.500000,
       132120.500000,
       66251.000000,
       265657.000000,
       133287.500000,
       546102.000000,
       274982.500000,
       1083604.000000,
       540816.500000,
       2160562.000000,
       1074971.000000,
       4423170.000000,
       2187084.000000,
       8907054.500000,
       4625695.000000,
       17826679.000000,
       9397566.000000,
       35653590.500000,
       18984357.500000,
       71369124.500000,
       37915319.000000],dtype="float32")



values_normal = values[::2]
values_optimized = values[1::2]

print(t)


def w_n(x) -> np.float32:
       return np.float(5*t[i])

for i in range(len(values_normal)):
       values_normal[i] = w_n(i) / values_normal[i]
       values_optimized[i] = w_n(i) / values_optimized[i]

# Show axes
fig, ax = plt.subplots()
ax.plot(t, values_normal)
ax.plot(t, values_optimized)

# Add textes
plt.text(0.45, 0.65, 'Vectorized version\n (Flags: -O3 -ffast-math -march=native)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
plt.text(0.50, 0.24, 'Unvectorized version (Flags: -O3 -fno-tree-vectorize)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
plt.text(0.30, 0.10, 'Unoptimized version (Flags: -O0)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
# Some settings
plt.grid()
ax.grid()

plt.suptitle('{} [{}]'.format(title,unity),fontsize=16,horizontalalignment='left', verticalalignment='bottom',x=0.125, y=.96)
plt.title("\nIntel, 11th Gen Intel(R) Core(TM) i9-11900H @ 2.50GH\n" +
              "L1d cache: 384 KiB\n" +
              "L1i cache: 256 KiB\n" +
              "L2 cache: 10 MiB\n" + 
              "L3 cache: 24 MiB\n" + 
              "Compiler : GCC 11.3.0  Flags : -O3 -mfma -fno-tree-vectorize\n", fontsize=8, loc='left')


ax.set_xlabel("Input size") 

plt.xscale("log")
 



ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

plt.show()