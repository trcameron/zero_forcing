#####################################################
#   Restricted Numerical Range Make Include File    #
#####################################################
SYSTEM = x86-64_osx
LIBFORMAT = static_pic
CPLEXDIR = /Applications/CPLEX_Studio201/cplex
CONCERTDIR = /Applications/CPLEX_Studio201/concert
NAUTY = /usr/local/Cellar/nauty/27r1/bin
PYTHON = python3.7
N = 3
NUM = 10
# ---------------------------------------------------------------------
# Compiler selection 
# --------------------------------------------------------------------- 
CCC = g++ -std=gnu++2a
# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCOPT = -Ofast -march=native -mtune=native
# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CPLEXBINDIR = $(CPLEXDIR)/bin/$(BINDIST)

CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) 
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR = $(CPLEXDIR)/include
LOCALINCDIR = $(PWD)/cpp/headers

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(LOCALINCDIR)