#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.lines  as mlines
import numpy as np
import csv
import math
import os
from scipy.interpolate import interp1d

size = 10000000 #20000
bondLength = 14

t=[]
x=[]
y=[]
z=[]
rid=[]
xOut=[]
yOut=[]
zOut=[]
ridOut=[]
xMF=[]
yMF=[]
zMF=[]
resIDMF=[]
cellNumberXMF=[]
cellNumberYMF=[]
cellNumberZMF=[]
betaMF=[]
bondLengths=[]
numberOfCGAtoms = 0
#cellVectorX=[39.970,0,0]
#cellVectorY=[-7.2383,25.9598,0]
#cellVectorZ=[-54.2489,-5.7914,675.7011]
cellVectorX=[39.9620, 0.0000,  0.7852]
cellVectorY=[-7.2369,25.9598, -0.1422]
cellVectorZ=[ 0.0000, 0.0000,677.9000]
aFCC=16.322
cellsX=30
cellsY=40
cellsZ=5
radius=6700
overlapTOgap = 0.8

mineralToCollagenRadiusRatio = 1
overlapMineralisationLengthRatio = 1  # 0 minerals only in the gap region, 1 minerals in the whole gap and the overlap length.

with open('align.tcl','w') as file:
        file.write("lappend auto_path ~/Downloads/la1.0\n")
        file.write("lappend auto_path ~/Downloads/orient\n")
        file.write("mol load pdb 3hr2.pdb\n")
        file.write("  package require Orient\n")
        file.write("namespace import Orient::orient\n\n")
        file.write("set sel [atomselect top all]\n")
        file.write("set I [draw principalaxes $sel]\n")
        file.write("set A [orient $sel [lindex $I 2] {0 0 1}]\n")
        file.write("$sel move $A\n")
        file.write("set I [draw principalaxes $sel]\n")
        file.write("set A [orient $sel [lindex $I 1] {0 1 0}]\n")
        file.write("$sel move $A\n")
        file.write("set I [draw principalaxes $sel]\n")
        file.write("set sel [atomselect top all]\n")
        file.write("$sel writepdb  0-3hr2-aligned.pdb\n")
        file.write("exit\n")

os.system("vmd -dispdev text -eofexit -e align.tcl")
os.system("vmd -dispdev text -eofexit -e InitialPoints.tcl")

with open('0-3hr2-AA.csv', 'r') as csvfile:
	plots = csv.reader(csvfile, delimiter=',')
	for row in plots:
		t.append(float(row[0]))
		rid.append(float(row[0]))
		x.append(float(row[1]))
		y.append(float(row[2]))
		z.append(float(row[3]))


for i in range(len(t)):
	for j in range(i+1,len(t)):
		if (z[i] > z[j]):
			ridDummy = rid[i]
			xDummy = x[i]
			yDummy = y[i]
			zDummy = z[i]
			rid[i] = rid[j]
			x[i] = x[j]
			y[i] = y[j]
			z[i] = z[j]
			rid[j] = ridDummy
			x[j] = xDummy
			y[j] = yDummy
			z[j] = zDummy

with open('0-3hr2-AA-Sorted.csv', 'w') as file:
	for i in range(len(t)):
		file.write("{:f},{:f},{:f},{:f}\n".format(t[i],x[i],y[i],z[i]))
			


t_new = np.linspace(min(t), max(t),size)

f1 = interp1d(t, x, kind='quadratic')
f2 = interp1d(t, y, kind='quadratic')
f3 = interp1d(t, z, kind='quadratic')
f4 = interp1d(t, rid, kind='quadratic')

x_new = f1(t_new)
y_new = f2(t_new)
z_new = f3(t_new)
rid_new = f4(t_new)

x1 = x_new[0]
y1 = y_new[0]
z1 = z_new[0]
rid1 = rid_new[0]


#xOut.append(x1)
#yOut.append(y1)
#zOut.append(z1)
outputLength = 0 #1
numberOfCGAtoms = 0 #1

flag = 0
for i in range(size):
	x2 = x_new[i]
	y2 = y_new[i]
	z2 = z_new[i]
	rid2 = rid_new[i]
	distance = math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
	if ((distance >= bondLength/2) and (flag == 0)):
		xOut.append(x2)
		yOut.append(y2)
		zOut.append(z2)
		ridOut.append(rid2)
		x1 = x2
		y1 = y2
		z1 = z2
		rid1 = rid2
		outputLength = outputLength + 1
		numberOfCGAtoms = numberOfCGAtoms + 1
		flag = 1
	if (distance >= bondLength):
		xOut.append(x2)
		yOut.append(y2)
		zOut.append(z2)
		ridOut.append(rid2)
		x1 = x2
		y1 = y2
		z1 = z2
		rid1 = rid2
		outputLength = outputLength + 1
		bondLengths.append(distance)
		numberOfCGAtoms = numberOfCGAtoms + 1

with open('1-mapping-information.txt', 'w') as file:
	file.write("{}\n".format(numberOfCGAtoms))
	file.write(" \n")
	file.write("[1,{}]\n".format(math.floor((ridOut[0]+ridOut[1])/2)))
	for i in range(1,outputLength-1):
		file.write("[{},{}]\n".format(math.ceil((ridOut[i]+ridOut[i-1])/2),math.floor((ridOut[i]+ridOut[i+1])/2)))
	file.write("[{},1054]\n".format(math.ceil((ridOut[outputLength-2]+ridOut[outputLength-1])/2)))

with open('1-3hr2-CG.pdb', 'w') as file:
	file.write("HEADER  COARSE GRAINED STRUCTURE FOR 3HR2\n")
	file.write("CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2\n")
	for i in range(outputLength):
		file.write("ATOM{:7d}  CA  CG1 A{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  1.00           C\n".format(i+1,1,xOut[i],yOut[i],zOut[i]))

with open('1-BondLengths.csv', 'w') as file:
	writer = csv.writer(file)
	for i in range(outputLength-1):
		row = [i,bondLengths[i]]
		writer.writerow(row)



with open('2-Microfibril-Square.pdb','w') as file:
	file.write("HEADER  CG model for collagen microfibril from 3HR2\n")
	atomID = 0
	cellNumberX = 0
	cellNumberY = 0
	cellNumberZ = 0
	resID = 0
	for i in range(cellsX):
		for j in range(cellsY):
			for k in range(cellsZ):
				cellNumberX = i + 1
				cellNumberY = j + 1
				cellNumberZ = k + 1
				resID = resID + 1
				for l in range(outputLength):
					x = xOut[l]+(-cellsX/2+i)*cellVectorX[0]+(-cellsY/2+j)*cellVectorY[0]+(-4+k)*cellVectorZ[0]
					y = yOut[l]+(-cellsX/2+i)*cellVectorX[1]+(-cellsY/2+j)*cellVectorY[1]+(-4+k)*cellVectorZ[1]
					z = zOut[l]+(-cellsX/2+i)*cellVectorX[2]+(-cellsY/2+j)*cellVectorY[2]+(-4+k)*cellVectorZ[2]
					if (z >= 0 and z <= 5*cellVectorZ[2]):
						atomID = atomID + 1
						file.write("ATOM{:7d}  CA  CG1 A{:4d}    {:8.3f}{:8.3f}{:8.3f} {:5.2f} {:5.2f}           C\n".format(atomID,resID,x,y,z,cellNumberX,cellNumberY))


with open('2-Microfibril-Square.gro','w') as file:
	file.write("CG model for collagen microfibril from 3HR2\n")
	file.write("{}\n".format(atomID))
	atomID = 0
	cellNumberX = 0
	cellNumberY = 0
	cellNumberZ = 0
	resID = 0
	for i in range(cellsX):
		for j in range(cellsY):
			for k in range(cellsZ):
				cellNumberX = i + 1
				cellNumberY = j + 1
				cellNumberZ = k
				resID = resID + 1
				betaMF.append(0)
				for l in range(outputLength):
					x = xOut[l]+(-cellsX/2+i)*cellVectorX[0]+(-cellsY/2+j)*cellVectorY[0]+k*cellVectorZ[0]
					y = yOut[l]+(-cellsX/2+i)*cellVectorX[1]+(-cellsY/2+j)*cellVectorY[1]+k*cellVectorZ[1]
					z = zOut[l]+(-cellsX/2+i)*cellVectorX[2]+(-cellsY/2+j)*cellVectorY[2]+k*cellVectorZ[2]
					if (z > 5*cellVectorZ[2]):
						#x = x - 5*cellVectorZ[0]
						#y = y - 5*cellVectorZ[1]
						z = z - 5*cellVectorZ[2]
					atomID = atomID + 1
					xMF.append(x/10)
					yMF.append(y/10)
					zMF.append(z/10)
					cellNumberXMF.append(cellNumberX)
					cellNumberYMF.append(cellNumberY)
					cellNumberZMF.append(cellNumberZ)
					resIDMF.append(resID)
					if ((x+50)*(x+50)+(y-150)*(y-150)<radius):
						betaMF[resID-1] = 1
					if (k == 0):
						file.write("{:5d}O{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 1):
						file.write("{:5d}A{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 2):
						file.write("{:5d}B{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 3):
						file.write("{:5d}C{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 4):
						file.write("{:5d}D{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 5):
						file.write("{:5d}E{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 6):
						file.write("{:5d}F{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 7):
						file.write("{:5d}G{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 8):
						file.write("{:5d}H{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 9):
						file.write("{:5d}I{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
					elif (k == 10):
						file.write("{:5d}J{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x/10,y/10,z/10))   # nm
	file.write("  0.00000   0.00000  337.85055   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000\n")

atoms = atomID
newAtoms = 0
for i in range(atoms):
	index = resIDMF[i]
	if (betaMF[index-1] == 1):
		newAtoms = newAtoms + 1


numberOfCollagenBeads = newAtoms
with open('3-Microfibril-Circular.gro','w') as file:
	file.write("CG for collagen microfibrill written by Mahdi Tavakol (mahditavakol90@gmail.com)\n")
	file.write("{}\n".format(numberOfCollagenBeads))
	atoms = atomID
	atomID = 0
	resID = 0
	currResID = 0
	for i in range(atoms):
		index = resIDMF[i]
		if (betaMF[index-1] == 1):
			k = cellNumberZMF[i]
			#resID = resIDMF[i]
			if (currResID != resIDMF[i]):
				currResID = resIDMF[i]
				resID = resID + 1
			cellNumberX = cellNumberXMF[i]
			cellNumberY = cellNumberYMF[i]
			atomID = atomID + 1
			x = xMF[i]
			y = yMF[i]
			z = zMF[i]
			if (k == 0):
				file.write("{:5d}O{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 1):
				file.write("{:5d}A{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 2):
				file.write("{:5d}B{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 3):
				file.write("{:5d}C{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 4):
				file.write("{:5d}D{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 5):
				file.write("{:5d}E{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 6):
				file.write("{:5d}F{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 7):
				file.write("{:5d}G{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 8):
				file.write("{:5d}H{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 9):
				file.write("{:5d}I{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
			elif (k == 10):
				file.write("{:5d}J{:02d}{:02d}   CL{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(resID,cellNumberX,cellNumberY,atomID,x,y,z))   # nm
	file.write("  0.00000   0.00000  337.85055\n")



# Mineralization
lastCollagenResID = resID
numberOfMineralBeads = 0

resID = resID + 1


with open('MinerTemp.gro','w') as file:
	R = 100.0
	nX = int(2*radius/aFCC)
	nY = int(2*radius/aFCC)
	nZ = int(5*cellVectorZ[2]/aFCC)
	basis1 = [ 0 , 0 , 0 ]
	basis2 = [ 0 ,0.5,0.5]
	basis3 = [0.5, 0 ,0.5]
	basis4 = [0.5,0.5, 0 ]
	basis = [basis1,basis2,basis3,basis4]
	for i in range(nX):
		for j in range(nY):
			for k in range(nZ):
				x = -50 - radius + i*aFCC
				y = 150 - radius + j*aFCC
				z = 0 + k*aFCC
				# resID = resID + 1
				cellZIndex = z%cellVectorZ[2]
				overlapSize = (overlapTOgap*cellVectorZ[2]) / (1+overlapTOgap)

				for base in basis:
					xFCC = x + aFCC*base[0]
					yFCC = y + aFCC*base[1]
					zFCC = z + aFCC*base[2]
					if ((x+50)*(x+50)+(y-150)*(y-150)<mineralToCollagenRadiusRatio * radius and cellZIndex >= (1-overlapMineralisationLengthRatio)*overlapSize):
						atomID = atomID + 1
						numberOfMineralBeads = numberOfMineralBeads + 1
						file.write("{:5d}Z{:04d}   HA11111{:8.3f}{:8.3f}{:8.3f}\n".format(resID,i,xFCC/10,yFCC/10,zFCC/10))   # nm



numberOfBeads = numberOfCollagenBeads + numberOfMineralBeads


with open('4-MT-Mineralized.gro','w') as file:
	file.write("CG for mineralized collagen microfibrill written by Mahdi Tavakol (mahditavakol90@gmail.com)\n")
	file.write("{}\n".format(numberOfBeads))
	MTFile = open('3-Microfibril-Circular.gro', 'r')
	next(MTFile)
	next(MTFile)
	beads = 0
	for line in MTFile:
		if ( beads < numberOfCollagenBeads):
			file.write(line)
			beads = beads + 1
	MinerFile = open('MinerTemp.gro', 'r')
	for line in MinerFile:
		file.write(line)
	file.write("  0.00000   0.00000  337.85055\n")
	MTFile.close()
	MinerFile.close()


cutoffs = [12.37,12.38,12.39,12.45,12.46]

for i in cutoffs:
	with open('badcontact.tcl','w') as file:
		distance = i
		file.write("mol load gro 4-MT-Mineralized.gro\n")
       		file.write("set sel [atomselect top all]\n")
        	file.write("$sel set beta 0\n")
        	file.write("set bad [atomselect top {resid > ")
        	file.write("{}".format(lastCollagenResID))
        	file.write(" and ( within ")
        	file.write("{}".format(distance))
        	file.write(" of resid < ")
        	file.write("{}".format(lastCollagenResID+1))
        	file.write(")}]\n")
        	file.write("$bad set beta 1\n")
        	file.write("set out [atomselect top {beta 0}]\n")
        	file.write("$out writepdb 4-MT-Mineralized-{}.pdb\n".format(distance))
		file.write("mol delete top\n")
		file.write("mol load pdb 4-MT-Mineralized-{}.pdb\n".format(distance))
       		file.write("set sel [atomselect top all]\n")
		file.write("$sel writegro 4-MT-Mineralized-{}.gro\n".format(distance))
		file.write("mol delete top\n")
		file.write("exit\n")
		file.close()
		os.system("vmd -dispdev text -eofexit -e badcontact.tcl")

