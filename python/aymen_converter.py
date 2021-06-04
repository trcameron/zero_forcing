# Aymen Converter
from sys import stdin
from nauty_geng_reader import graph6
from networkx import Graph
import traceback

###############################################
###             main                        ###
###############################################
def main():
    try:
        # read input stream
        for line in stdin:
            # build graph
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            # write graph
            f = open("graphs%d.txt"%g.order(),"a+")
            f.write("%d\n"%g.order())
            for e in g.edges:
                f.write("%d %d\n"%(e[0],e[1]))
            f.write("-1\n")
            f.close()
    except Exception as e:
        print(traceback.format_exc())
if __name__ == '__main__':
    main()