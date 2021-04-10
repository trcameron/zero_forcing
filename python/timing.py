# Timing
from matplotlib import pyplot as plt

wave_time = [3.704E-4,9.960E-3,7.779E-2,6.411E-1]
ga_time = [5.284E-2,3.647E-1,7.937E-1,1.589E+0]
dom = [10+10*k for k in range(4)]

plt.semilogy(dom,wave_time)
plt.semilogy(dom,ga_time)
plt.savefig("timing.png",dpi=400)