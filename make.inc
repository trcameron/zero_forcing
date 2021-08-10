#####################################################
#   Zero-Forcing Make Include File    				#
#####################################################
BOOSTDIR = /usr/local/Cellar/boost/1.76.0
GUROBIDIR = /Library/gurobi912/mac64
NAUTY = /usr/local/Cellar/nauty/27r1/bin
ORDER = 3
# ---------------------------------------------------------------------
# Compiler selection 
# --------------------------------------------------------------------- 
CCC = g++ -std=c++11
# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCOPT =
# ---------------------------------------------------------------------
# Include directories
# ---------------------------------------------------------------------
BOOSTINCDIR = $(BOOSTDIR)/include
GUROBILIBDIR = $(GUROBIDIR)/lib
GUROBIINCDIR = $(GUROBIDIR)/include
LOCALINCDIR = $(PWD)/cpp/headers
# ---------------------------------------------------------------------
# Link and include flags
# ---------------------------------------------------------------------
CCLNFLAGS = -L$(GUROBILIBDIR) -lgurobi91 -lgurobi_c++
CCINCFLAGS = -I$(BOOSTINCDIR) -I$(GUROBIINCDIR) -I$(LOCALINCDIR)