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

filename = "EE_3D_%sx_%sy_%sz_L%d_mass%s.txt" %(bc_x,bc_y,bc_z,L,decimalStr(mass))
fout = open(filename, 'w')

alpha = [1, 2]
LA_x_max = (L)/2
LA_y     = L
LA_z     = L

#Calculate all needed correlators:
X_from0, P_from0 = free_boson_3D.getXP_from0(L, LA_x_max, LA_y, LA_z, bc_x, bc_y, bc_z, mass)

#Loop over all cylinders:
for LA_x in range(1,LA_x_max+1):
  t1_A = time.clock() #for timing
  print "\nLA = %d" %LA_x
  sys.stdout.flush()
  XA,PA = free_boson_3D.getXAPA_fromXP0(L, LA_x, LA_y, LA_z, X_from0, P_from0, mass)
  
  S_alpha_Evs = free_boson_3D.getEntropy_Evs(XA,PA)
  S_alpha     = free_boson_3D.getEntropy(S_alpha_Evs,alpha)
  
  #Save results to file:
  fout.write("%d" %LA_x)
  for Sn in S_alpha:
    fout.write(" %.15f" %Sn)
  fout.write("\n")
  fout.flush()

  t2_A = time.clock()
  print "\nTime for LA = %d: %.5f sec. " %(LA_x, t2_A-t1_A)
  sys.stdout.flush()
#End loop over cylinders

t2 = time.clock()
print "Total elapsed time: %.5f sec." %(t2-t1)
