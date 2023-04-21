import numpy as np 
import matplotlib.pyplot as plt

plt.style.use('bmh')

title = "Test plot"
flags = "-03"
unity = "flops / cycles"
output_path = "sample_plot.png"

if unity == "" or flags == "" or title == "" or output_path == "":
       print("Missing something")
       exit()

# Plot values
scalar = [0,3,3,3]
vector = [0,3,12,12]
ref = [0,0.125,0.5,32]


# Show axes
fig, ax = plt.subplots()
ax.plot(ref, scalar,label='Scalar roofline')
ax.plot(ref, vector,label='Vector roofline')
ax.plot(0.291, 2.8, 'go', label='Comparaison 1 - b') 
ax.plot(0.291, 2, 'ro', label='Comparaison 2 - b') 
ax.plot(0.291, 1.17, 'bo', label='Comparaison 3 - b') 
ax.plot(0.076, 1.83, 'yo', label='Comparaison 1/2 - d')        
ax.plot(0.076, 1.17, 'co', label='Comparaison 3 - d') 


# Some settings
plt.grid()
ax.grid()
plt.title('Roofline plot')
ax.set_xlabel('Intensity (W/Q) [flops/bytes]')
ax.set_ylabel('Performance (W/T) [flops/cycle]')

ax.set_xlim(left=1/32, right=32)
ax.set_ylim(bottom = 1/32, top=32)

plt.legend()

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xscale("log", base=2);
ax.set_yscale("log", base=2);

ticks= [pow(2,i) for i in range(-5,5)]


plt.xticks(ticks=ticks)
plt.yticks(ticks=ticks)

# Save figure
fig.savefig(output_path)
plt.show()
