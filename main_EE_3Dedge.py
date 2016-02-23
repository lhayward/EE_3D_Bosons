import argparse
import numpy as np
import os.path  #to check if file exists
import sys  #for sys.stdout.flush()
import time
import matplotlib.pyplot as plt

import free_boson_3D

########################################  decimalStr  #########################################
# Converts a number to a string and replaces all commas with hyphens
# (used when naming files, where we don't want periods in the file name)
###############################################################################################
def decimalStr(num):
  res = str(num)
  length = len(res)
  index = res.find('.')
  if index >= 0:
    res = res[0:index] + '-' + res[(index+1):length]
  return res
# end decimalStr(num) function

#########################################  readArray  #########################################
# Takes a string representation of an array and removes the '[' and ']' characters
###############################################################################################
# def readArray(line):
#   start = max( 0, line.find('[') )
#   end = min( len(line), line.find(']') )
#   return line[start+1:end] 
####### end readArray(line) function #######

########################################  readParams  #########################################
# Reads in values from input file. Input should have the following form:
#
# L    = ___ (int)
# bc_x = ___ (string: 'PBC' or 'APBC')
# bc_y = ___ (string: 'PBC' or 'APBC')
# bc_z = ___ (string: 'PBC' or 'APBC')
# mass = ___ (float)
###############################################################################################
def readParams(filename):
  L    = 1
  bc_x = 'PBC'
  bc_y = 'PBC'
  bc_z = 'PBC'
  massterm = 1
  if os.path.isfile(filename):
    fin = open(filename,'r')
    
    line=fin.readline()
    L = int(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    bc_x = line[ max(0,line.find('='))+1:].strip()
    
    line=fin.readline()
    bc_y = line[ max(0,line.find('='))+1:].strip()
    
    line=fin.readline()
    bc_z = line[ max(0,line.find('='))+1:].strip()
    
    line=fin.readline()
    mass = float(line[ max(0,line.find('='))+1:])
    
    fin.close()
  return L, bc_x, bc_y, bc_z, mass

###############################################################################################
###########################################  main  ############################################
###############################################################################################

parser=argparse.ArgumentParser(description="Code to calculate EE for 3D free bosons")
parser.add_argument('-f', '--file')
args=parser.parse_args()

###### Read input from file: ######
inFile = "input"
if args.file != None:
  inFile = inFile + "_" + args.file
inFile = inFile + ".txt"
print "Input file: %s" %inFile

L, bc_x, bc_y, bc_z, mass = readParams(inFile)
###################################

print "L          = %d" %L
print "BC along x = %s" %bc_x
print "BC along y = %s" %bc_y
print "BC along z = %s" %bc_z
print "mass       = %f" %mass
sys.stdout.flush()

t1 = time.clock() #for timing

filename = "EE_evals_3D_%sx_%sy_%sz_L%d_mass%s.txt" %(bc_x,bc_y,bc_z,L,decimalStr(mass))
fout = open(filename, 'w')

LA_x = L/2
LA_y = L/2
LA_z = L

#XA,PA = free_boson_3D.getXAPA_edge_K(L,LA_x,LA_y,LA_z,bc_x,bc_y,bc_z,mass)
XA,PA = free_boson_3D.getXAPA_edge_corr(L,LA_x,LA_y,LA_z,bc_x,bc_y,bc_z,mass)

t2 = time.clock()
print "\nTime to build XA, PA: " + str(t2-t1) + " sec."
sys.stdout.flush()

S_alpha_Evs = free_boson_3D.getEntropy_Evs(XA,PA)

#Save results to file:
for ev in S_alpha_Evs:
  fout.write(" %.20f\n" %ev)
#fout.flush()
fout.close()

t2 = time.clock()
print "Total elapsed time: " + str(t2-t1) + " sec."
