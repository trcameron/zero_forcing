# Zero Forcing Test
from sys import stdin
from time import time
from networkx import Graph
from nauty_geng_reader import graph6
from matplotlib import pyplot as plt
from Networkx import zero_forcing, psd_rule
from zero_forcing_ip import zf_std

###############################################
###             main                        ###
###############################################
def main():
    try:
        # testing parameters
        et_ga = []
        et_ip = []
        # read input stream
        for line in stdin:
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            
            start_time = time()
            zf1 = zero_forcing(g)
            end_time = time()
            et_ga.append(end_time-start_time)
            
            start_time = time()
            zf2,_,_,_ = zf_std(a)
            end_time = time()
            et_ip.append(end_time - start_time)
        # plot testing results
        plt.subplot(111)
        plt.plot(et_ga)
        plt.plot(et_ip)
        plt.legend(('ga_time','ip_time'))
        plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()