# Zero Forcing Test
from sys import stdin
from time import time
from networkx import Graph, draw_spectral
from nauty_geng_reader import graph6
from matplotlib import pyplot as plt
from Networkx import zero_forcing as zf_gen
from zero_forcing_ip import zero_forcing as zf_ip

###############################################
###             main                        ###
###############################################
def main():
    try:
        # read input stream
        for line in stdin:
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            draw_spectral(g,ax=plt.subplot(121))
            
            start_time = time()
            zf1 = zf_gen(g)
            end_time = time()
            et1 = (end_time-start_time)
            
            start_time = time()
            zf2,_,_,_ = zf_ip(a)
            end_time = time()
            et2 = (end_time - start_time)
            
            plt.subplot(122)
            plt.axis("off")
            plt.text(0.5,1.0,"zf_gen = %d, time = %.4e"%(zf1,et1),size=12,ha="center")
            plt.text(0.5,0.5,"zf_ip = %d, time = %.4e"%(zf2,et2),size=12,ha="center")
            plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()