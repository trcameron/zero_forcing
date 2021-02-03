# Tree Paramters
from nauty_geng_reader import graph6, sparse6
from networkx import Graph, draw, all_simple_paths, connected_components
from matplotlib import pyplot as plt
from sys import stdin

###############################################
###             nylen min rank              ###
###############################################
def nylen_min_rank(g):
    count = 0
    while(len(g.nodes)>0):
        # degree of each vertex
        deg = g.degree
        # leafs in graph
        leaf = []
        for k in g.nodes:
            if(deg[k]==1):
                leaf.append(k)
        # Nylen path
        nylen = False
        while(not(nylen)):
            target = leaf.pop()
            source = leaf.pop()
            for x in all_simple_paths(g,source,target):
                x_deg = g.degree(x)
                high_deg_ind = [k for k in range(len(x)) if x_deg[k]>=3]
                if(len(high_deg_ind)<=1):
                    nylen = True
                    break
        # update count
        if(len(high_deg_ind)==1):
            count += 2
            g.remove_node(high_deg_ind[0])
        else:
            count += len(x) - 1
            g.remove_nodes_from(x)
    # return minimum rank
    return count
###############################################
###             main                        ###
###############################################
def main():
    try:
        a = sparse6(bytes(":Ccf",'utf-8'))
        g = Graph(a)
        nylen_min_rank(g)
        draw(g,ax=plt.subplot(111))
        plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()
