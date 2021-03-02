# Zero Forcing Test
from os import popen
from time import time
from networkx import Graph, draw
from nauty_geng_reader import graph6
from matplotlib import pyplot as plt
from Networkx import zero_forcing, psd_rule
from zero_forcing_ip import zf_std

start_ord = 5
itnum = 5
num = 100
###############################################
###             main                        ###
###############################################
def main():
    try:
        et_ga_avg = []
        et_ip_avg = []
        fail_count = 0
        n = start_ord
        for k in range(itnum):
            print('/usr/local/Cellar/nauty/27r1/bin/genrang -g %d %d'%(n,num))
            f = popen('/usr/local/Cellar/nauty/27r1/bin/genrang -g %d %d'%(n,num))
            line = f.readline()
            
            et_ga = []
            et_ip = []
            while line != '':
                a = graph6(bytes(line.rstrip(),'utf-8'))
                g = Graph(a)
                
                start_time = time()
                zf1 = zero_forcing(g)
                end_time = time()
                et_ga.append(end_time-start_time)
                
                start_time = time()
                zf2 = zf_std(a)
                end_time = time()
                et_ip.append(end_time - start_time)
                
                if(zf1[0]!=round(zf2[0])):
                    color_map = []
                    fail_count += 1
                    for i in range(len(zf2[1])):
                        if(round(zf2[1][i])==1):
                            color_map.append('#0000FF')
                        else:
                            color_map.append('#C0C0C0')
                    draw(g,with_labels=True,node_color=color_map,ax=plt.subplot(121))
                    plt.subplot(122)
                    plt.axis("off")
                    plt.text(0.5,0.75,"GA ZF-Number = %d"%zf1[0],size=12,ha="center")
                    plt.text(0.5,0.25,"IP ZF-Number = %.2f"%zf2[0],size=12,ha="center")
                    plt.savefig("rand_test/ga_fail%d.png"%fail_count,dpi=400)
                    plt.cla(); plt.clf(); plt.close()
                    
                line = f.readline()
            
            n = n*2
                
            et_ga_avg.append(sum(et_ga)/len(et_ga))
            et_ip_avg.append(sum(et_ip)/len(et_ip))
            
        plt.subplot(111)
        dom = [start_ord*2**k for k in range(itnum)]
        plt.semilogy(dom,et_ga_avg)
        plt.semilogy(dom,et_ip_avg)
        plt.legend(('ga_avg_time','ip_avg_time'))
        plt.savefig("rand_test/avg_time.png",dpi=400)
        plt.cla(); plt.clf(); plt.close()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()