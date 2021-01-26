# Nauty GENG Reader
from sys import stdin
from numpy import zeros
from networkx import Graph, draw
from matplotlib import pyplot as plt

def graph6(bytes_in):
    def data_to_n(data):
        """Read initial one-, four- or eight-unit value from digraph6 integer sequence.
        Return (value, rest of seq.)"""
        if data[0] <= 62:
            return data[0], data[1:]
        elif data[1] <= 62:
            return (data[1] << 12) + (data[2] << 6) + data[3], data[4:]
        else:
            return (data[2] << 30) + (data[3] << 24) + (data[4] << 18) + (data[5] << 12) + (data[6] << 6) + data[7], data[8:]
    def bits():
        """Returns sequence of individual bits from 6-bit-per-value list of data values."""
        for d in data:
            for i in [5, 4, 3, 2, 1, 0]:
                yield (d >> i) & 1
    # subtract 63 from each bit, check for bits that are over 126
    data = [c - 63 for c in bytes_in]
    if any(c > 63 for c in data):
        raise Exception("graph6 characters must be in range(63, 127)")
    # extract size of graph (n) and remaining bits which hold edge information
    n, data = data_to_n(data)
    # check if we have the correct number of bits
    nd = (n*(n-1)//2 + 5) // 6
    if len(data) != nd:
        raise Exception(f"Expected {nd*6} bits and got {len(data) * 6} in graph6")
    # build adjacency matrix    
    a = zeros((n,n),dtype=float)
    for (i, j), b in zip([(i, j) for j in range(1,n) for i in range(j)], bits()):
        if b:
            a[i,j] = 1; a[j,i] = 1
    # return
    return a
###############################################
###             main                        ###
###############################################
def main():
    try:
        # read input stream
        for line in stdin:
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            draw(g,ax=plt.subplot(111))
            plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()