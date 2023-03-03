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
t = [0,1,2,3,4,5]
s = [0,1,2,3,4,5]
s2 = [0,0.5,1.2,2,2,2.3]

# Show axes
fig, ax = plt.subplots()
ax.plot(t, s)
ax.plot(t, s2)

# Add textes
plt.text(0.55, 0.63, 'text', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
plt.text(0.65, 0.42, 'text2', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
# Some settings
plt.grid()
ax.grid()
ax.xaxis.set_label_position('top')
plt.title('{} [{}]'.format(title,unity),fontsize=16, loc='left')
ax.set_xlabel("Flags : {}".format(flags), loc='left')
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

# Save figure
fig.savefig(output_path)
plt.show()