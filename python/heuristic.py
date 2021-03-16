# Heuristic Algorithms for the Zero-Forcing Number
from sys import stdin
from nauty_geng_reader import graph6
from networkx import Graph, draw
from matplotlib import pyplot as plt
from numpy import array, sum, zeros
from collections import deque
import traceback

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
                colored[v] = 1
                # update the number of colored neighbors for each vertex in s that is a neighbor of v. 
                for w in g.neighbors(v):
                    if(colored[w]):
                        count[w] += 1
                        if(count[w]==(g.degree[w]-1)):  # update vertices that can do forcing
                            stack.append(w)
                # update the number of colored neighbors of v
                count[v] = sum([colored[u] for u in g.neighbors(v)])
                if(count[v]==(g.degree[v]-1)):
                    stack.append(v)
    # return
    return set([v for v in g.nodes if colored[v]])
###############################################
###             add                         ###
###############################################
def add(v,c):
    return set([v])
###############################################
###             f                           ###
###############################################
def f(q,v,c):
    return len(q)
###############################################
###             heuristic                   ###
###############################################
def heuristic(g):
    # initialize colored set and closure
    z = set()
    c = set()
    # while z is not a zero-forcing set
    while(c != set(g.nodes)):
        c_new = set()
        a = set()
        for v in (set(g.nodes)-c):
            s = c.union(add(v,c))
            if(f(closure(g,s),v,c)>f(c_new,v,c)):
                a = add(v,c)
                c_new = closure(g,s)
        c = c_new.copy()
        z = z.union(a)
        for v in z:
            if(closure(g,z-set([v]))==g.nodes):
                z = z-set([v])
    # return
    return z
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
            print(z)
            draw(g,with_labels=True,ax=plt.subplot(111))
            plt.show()
    except Exception as e:
        print(traceback.format_exc())
if __name__ == '__main__':
    main()