# Wavefront Algorithm for the PSD Zero-Forcing Number
from sys import stdin
from nauty_geng_reader import graph6
from networkx import Graph, draw, path_graph
from networkx.algorithms.simple_paths import all_simple_paths
from networkx.classes.function import is_frozen
from matplotlib import pyplot as plt
from numpy import array, sum, zeros
from collections import deque

###############################################
###             closure                     ###
###############################################
def closure(g,s):
    # graph order
    n = len(g.nodes)
    # initialize colored set, subgraph of white vertices, and active stack
    colored = s.copy()
    subg = g.subgraph(set(g.nodes()) - colored)
    active = deque()
    for u in colored:
        nbhd = set(g.neighbors(u))  # nbhd of u
        wnbhd = nbhd-colored        # white neighbors of u
        for w in wnbhd:
            paths = []
            for v in (wnbhd-set([w])):
                paths += all_simple_paths(subg,w,v)       # all white wv paths
            if(len(paths)==0):
                active.append([u,w])
    # while there are still active vertices
    while(active):
        a = active.pop()
        if(not (a[-1] in colored)):
            # force coloring
            colored.add(a[-1])
            # update subg
            subg = g.subgraph(set(g.nodes()) - colored)
            # update active stack
            for u in colored:
                nbhd = set(g.neighbors(u))  # nbhd of u
                wnbhd = nbhd-colored        # white neighbors of u
                for w in wnbhd:
                    paths = []
                    for v in (wnbhd-set([w])):
                        paths += all_simple_paths(subg,w,v)       # all white wv paths
                    if(len(paths)==0):
                        active.append([u,w])
    # return colored
    return colored
###############################################
###             main                        ###
###############################################
def main():
    g = path_graph(5)
    s = set([2])
    s_cl = closure(g,s)
    print(s_cl)
if __name__ == '__main__':
    main()