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
t = np.arange(100,1600,100)

no_opt = np.array([9102870,
       72952522,
       245319576,
       601896191,
       1181908811,
       2036199002,
       3307251007,
       4885165963,
       6998298818,
       9422353610,
       12570256347,
       16620476517,
       20947859176,
       27077661363,
       33689801895],dtype="float32")

opt = np.array([2821972,
       28379393,
       90447904,
       209945378,
       415748272,
       740861830,
       1295163388,
       1861470116,
       2649276680,
       3817256482,
       5224313878,
       6534126738,
       8439327276,
       11704734476,
       13343157364],dtype="float32")

vector = np.array([718446,
       5833809,
       21627788,
       56107521,
       111089481,
       193568093,
       311542591,
       496585326,
       712785793,
       964874806,
       1300715581,
       1761700038,
       2456036513,
       3189672380,
       4095521249],dtype="float32")


def w_n(x) -> np.float32:
       return np.float32(2*(100*x+100)**3 + (100*x+100)**2)

for i in range(len(no_opt)):
       no_opt[i] = w_n(i) / no_opt[i]
       opt[i] = w_n(i) / opt[i]
       vector[i] = w_n(i) / vector[i] 
# Show axes
fig, ax = plt.subplots()
ax.plot(t, no_opt)
ax.plot(t, opt)
ax.plot(t, vector)

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
              "Compiler : GCC 11.3.0 \n", fontsize=8, loc='left')


ax.set_xlabel("Input size") 



ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

plt.show()