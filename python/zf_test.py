# Zero-Force Testing
from sys import stdin
from heuristic import heuristic
from wavefront import wavefront
from zero_forcing_ip import zf_std
from nauty_geng_reader import graph6
from networkx import Graph, draw
from matplotlib import pyplot as plt
###############################################
###             main                        ###
###############################################
def main():
    try:
        # read input stream
        for line in stdin:
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            z = heuristic(g)
            zf1 = len(z)
            
            zf2 = wavefront(g)
            zf3,s,x,y = zf_std(g)
            if(zf2 != zf3):
                color_map = []
                for i in range(len(s)):
                    if(round(s[i])==1):
                        color_map.append('#0000FF')
                    else:
                        color_map.append('#C0C0C0')
                draw(g,with_labels=True,node_color=color_map,ax=plt.subplot(121))
                plt.subplot(122)
                plt.axis("off")
                plt.text(0.5,0.75,"Heuristic = %d"%zf1,size=12,ha="center")
                plt.text(0.5,0.5,"Wavefront = %d"%zf2,size=12,ha="center")
                plt.text(0.5,0.25,"IP = %d"%zf3,size=12,ha="center")
                plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()