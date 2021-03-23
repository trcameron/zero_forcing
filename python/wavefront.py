# Wavefront Algorithm for the Zero-Forcing Number
from sys import stdin
from nauty_geng_reader import graph6
from networkx import Graph, draw
from matplotlib import pyplot as plt
from numpy import array, sum, zeros
from collections import deque

###############################################
###             closure                     ###
###############################################
def closure(g,s):
    # graph order
    n = len(g.nodes)
    # initialize colored and count arrays and stack
    colored = zeros(n,dtype=bool)
    count = zeros(n,dtype=int)
    stack = deque()
    for v in s:
        colored[v] = True  # only the vertices in s are colored
    for v in s:
        count[v] = sum([colored[u] for u in g.neighbors(v)])    # count the number of colored neighbors for each vertex in s
        if(count[v]==(g.degree[v]-1)):  # vertices that can do forcing
            stack.append(v)
    # while there are still vertices that can do forcing
    while(stack):
        u = stack.pop() # vertex that forces
        for v in g.neighbors(u):
            if(not colored[v]):  # vertex that is forced
                # color v
                colored[v] = True
                # update the number of colored neighbors for each colored vertex that is a neighbor of v. 
                for w in g.neighbors(v):
                    if(colored[w]):
                        count[w] += 1
                        if(count[w]==(g.degree[w]-1)):  # update vertices that can do forcing
                            stack.append(w)
                # update the number of colored neighbors of v
                count[v] = sum([colored[w] for w in g.neighbors(v)])
                if(count[v]==(g.degree[v]-1)):
                    stack.append(v)
                # break out of for v loop
                break
    # return
    return frozenset([v for v in g.nodes if colored[v]])
###############################################
###             wavefront                   ###
###############################################
def wavefront(g):
    # graph order
    n = len(g.nodes)
    # closure pairs
    c = set([(frozenset(),0)])
    # iterate over all possible cardinalities
    for k in range(1,n+1):
        # iterate over all closure pairs
        for (s,r) in c:
            # iterate over all vertices
            for v in g.nodes:
                nbhd = set(g.neighbors(v))
                sv = set([v])
                s_cl = closure(g,s.union(nbhd).union(sv))
                r_cl = r + len(sv-s) + max(len(nbhd-s)-1,0)
                if(r_cl<=k and not any([(s_cl,i) in c for i in range(k+1)])):
                    c = c.union(set([(s_cl,r_cl)]))
                    if(len(s_cl)==n):
                        return r_cl
###############################################
###             main                        ###
###############################################
def main():
    try:
        # read input stream
        for line in stdin:
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            zf = wavefront(g)
            print(zf)
            draw(g,with_labels=True,ax=plt.subplot(111))
            plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()