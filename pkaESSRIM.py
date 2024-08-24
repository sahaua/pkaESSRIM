## This program calculates and plots pka energy spectra from the output of SRIM-2013 code.
## Input file is the collision.txt file obtained from SRIM-2013.
## Recoil energy bins have to be provided in a file named as RecoilEnergyBins. 
## Author ---- Uttiyoarnab Saha ---- ## 
## Kolkata, 15.11.2021

##--------------------------------------------------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt

ifilename = input("Enter input file name: ")	#input("Enter input file name: ")
ofilename = input("Enter output file name: ")

ifile = open(ifilename, "r", encoding='latin1')
ifile1 = open('RecoilEnergyBins')
ofile = open(ofilename, "w")

def damageenergy(Z1,A1,Z2,A2,Er):
	twothd = 0.666666667
	threeq = 0.75
	sixth = 0.166666667
	onep5 = 1.5
	c1 = 30.724
	c2 = 0.0793
	c3 = 3.4008
	c4 = 0.40244
	dont = 40.0
	el = c1*Z1*Z2*(np.sqrt(Z1**twothd+Z2**twothd))*(A1+A2)/A2
	rel = 1/el
	fl1 = c2*Z1**twothd*(np.sqrt(Z2))*(A1+A2)**onep5
	fl2 = ((Z1**twothd+Z2**twothd)**threeq*A1**onep5*(np.sqrt(A2)))
	fl = fl1/fl2
	ep = Er*rel
	Tdam = Er/(1+fl*(c3*ep**sixth+c4*ep**threeq+ep))
	return(Tdam) 

def calcdefects(Tdam,Ed):
	numerator = 0.8*Tdam
	denominator = 2*Ed
	dpa = numerator/denominator
	return(dpa)

elem = []

ionNum = []
Eion = []
RecName = []
Erec = []

title = input('Enter Calculation Name: ')
print()
print('---- Incident Ion Details ----')
name_iion = input('Enter Incident Ion Symbol: ')
Nions = int(input('Total Number of Ions Simulated: '))
Einc = float(input('Enter Incident Ion Energy (keV): '))

print()
print('---- Target Details ----')
Ntargelm = int(input('Number of Elements in Target: '))
for i in range(Ntargelm):
	elem.append(input('Enter Element Symbol: '))

while True:
	line = ifile.readline()
	if (len(line.split(' ')) >= 3):
		if (line.split(' ')[2] == 'NOTES:'):
			for i in range(8):
				ifile.readline()
			break

while True:	
	line = ifile.readline()
	
	if line[0] == '=':
		for i in range(11):
			ifile.readline()
		line = ifile.readline()
		
	data = [(x) for x in line.split('³')]
	
	if len(data) == 13:
		ionNum.append(int(data[1]))
		Eion.append(float(data[2]))
		RecName.append(data[7].strip())
		Erec.append(float(data[8])*1E-03)
	
	if line == '':
		break	

ifile.close()

# Energy Bins for the Recoils  
n_recEn = int(ifile1.readline())
ifile1.readline()
recEnbin = np.zeros(n_recEn)

for i in range(n_recEn):
	recEnbin[i] = float(ifile1.readline())

ifile1.close()
	
countpka = np.zeros((n_recEn, Ntargelm))

for i in range(len(Erec)):
	for j in range(n_recEn-1):
		if (recEnbin[j] < Erec[i] and Erec[i] <= recEnbin[j+1]):
			jscore = j
			break
	for k in range(Ntargelm):
		if (elem[k] == RecName[i]):
			countpka[jscore][k] += 1
			break
	if (Erec[i] == recEnbin[0]):
		for k in range(Ntargelm):
				if (elem[k] == RecName[i]):
					countpka[0][k] += 1

countpka = countpka/Nions


countpkaNRT = np.zeros((n_recEn, Ntargelm)) 	# NRT Weighted PKA spectrum

Z1 = 28 #int(input('Enter Ion Z: '))
A1 = 58 #float(input('Enter Ion A: '))
Z2 = 28 #int(input('Enter Target Z: '))
A2 = 58 #float(input('Enter Target A: '))
Ed = 40 #float(input('Enter threshold lattice displacement energy: '))



for k in range(Ntargelm):
	s = 0
	s1 = 0
	for j in range(n_recEn-1):
		perkeV = recEnbin[j+1] - recEnbin[j]
		countpka[j][k] = countpka[j][k]/perkeV
		s = s + countpka[j][k]
		
		Tmid = 1E+03*0.5*(recEnbin[j] + recEnbin[j+1])	# in eV
		Tdmg = damageenergy(Z1,A1,Z2,A2,Tmid)
		Ndef = calcdefects(Tdmg,Ed)
		countpkaNRT[j][k] = countpka[j][k]*Ndef
		s1 = s1 + countpkaNRT[j][k]
		
	for j in range(n_recEn-1):
		countpka[j][k] = countpka[j][k]/s
		countpkaNRT[j][k] = countpkaNRT[j][k]/s1

# printing output

print('-------- PKA calculation for ', title, ' --------', file = ofile)

print('PKA Energy (keV)		PKA distribution', file = ofile)

for k in range(Ntargelm):
	print('		', elem[k], file = ofile)
	for j in range(n_recEn-1):
		print(recEnbin[j], '-', recEnbin[j+1], '	', '{:.6E}'.format(countpka[j][k]), file = ofile)
print('------------------------------------------------------------------', file = ofile)

dpivalue = 150
plt.figure(1, dpi = dpivalue)
plt.xscale('log')
plt.xlabel('$E_R$ (keV)', fontsize = 12)
plt.ylabel('PKA spectrum', fontsize = 12)
plt.step(recEnbin, countpka, where = 'post')
plt.title(title)
plt.legend(elem)
plt.show()

pka_LTEr = np.zeros((n_recEn,Ntargelm))
pka_LTErNRT = np.zeros((n_recEn,Ntargelm))

xdata = []
xdata = recEnbin[1:]

plt.figure(2, dpi = dpivalue)
plt.xscale('log')
plt.xlabel('$E_R$ (keV)', fontsize = 12)
plt.ylabel('Cumulative Fraction of PKA less than $E_R$', fontsize = 12)

print('Cumulative PKA fraction below E_R', file = ofile)

charc_Thalf = np.zeros(Ntargelm)
charc_ThalfNRT = np.zeros(Ntargelm)

for k in range(Ntargelm):
	flag = 0
	flag0 = 0
	ydata = []
	s = 0
	s0 = 0
	for j in range(n_recEn-1):
		s = s + countpka[j][k] * (recEnbin[j+1] - recEnbin[j])
		s0 = s0 + countpkaNRT[j][k] * (recEnbin[j+1] - recEnbin[j])
	s1 = 0
	s01 = 0
	for j in range(n_recEn-1):
		s1 = s1 + countpka[j][k] * (recEnbin[j+1] - recEnbin[j])
		s01 = s01 + countpkaNRT[j][k] * (recEnbin[j+1] - recEnbin[j])
		pka_LTEr[j+1][k] = s1/s
		pka_LTErNRT[j+1][k] = s01/s0
	
	print('		', elem[k], file = ofile)
	
	for j in range(n_recEn):
		if (pka_LTEr[j][k] >= 0.5 and flag == 0):
			charc_Thalf[k] = recEnbin[j]
			flag = 1
		
		if (pka_LTErNRT[j][k] >= 0.5 and flag0 == 0):
			charc_ThalfNRT[k] = recEnbin[j]
			flag0 = 1
		
		if j > 0:
			ydata.append(pka_LTEr[j][k])
		print(recEnbin[j], '{:.6E}'.format(pka_LTEr[j][k]), sep = ' ', file = ofile)

	plt.plot(xdata, ydata)	#,  where = 'pre'
	#plt.legend(elem[k])
plt.title(title)
plt.legend(elem)
plt.show()

# average energy of pka

pka_avgEn = np.zeros((Ntargelm))
for k in range(Ntargelm):
	s1 = 0
	s2 = 0
	for j in range(n_recEn-1):
		s1 = s1 + 0.5*(recEnbin[j] + recEnbin[j+1])*countpka[j][k] * (recEnbin[j+1] - recEnbin[j])
		s2 = s2 + countpka[j][k] * (recEnbin[j+1] - recEnbin[j])
	pka_avgEn[k] = s1/s2
print('-------------------------------------------------------------------', file = ofile)
print('Average energy of PKA:', file = ofile)
print(elem, pka_avgEn, sep = ' ', file = ofile)
print('', file = ofile)
print('Charactestic energy of PKA (E_R_1/2):', file = ofile)
print('Less than or equal to ', file = ofile)
print(elem, charc_Thalf, sep = ' ', file = ofile)
print('------------------------------------------------------------------', file = ofile)
print('>>>>>>=================================================<<<<<<', file = ofile)

print('PKA Energy (keV)		PKA distribution (NRT Weighted)', file = ofile)

for k in range(Ntargelm):
	print('		', elem[k], file = ofile)
	for j in range(n_recEn-1):
		print(recEnbin[j], '-', recEnbin[j+1], '	', '{:.6E}'.format(countpkaNRT[j][k]), file = ofile)
print('------------------------------------------------------------------', file = ofile)

print('', file = ofile)
print('Cumulative PKA fraction below E_R (NRT Weighted)', file = ofile)

for k in range(Ntargelm):
	print('		', elem[k], file = ofile)
	for j in range(n_recEn):
		print(recEnbin[j], '{:.6E}'.format(pka_LTErNRT[j][k]), sep = ' ', file = ofile)
		
print('', file = ofile)
print('Charactestic energy of PKA (E_R_1/2):', file = ofile)
print('Less than or equal to ', file = ofile)
print(elem, charc_ThalfNRT, sep = ' ', file = ofile)
print('------------------------------------------------------------------', file = ofile)

plt.figure(3, dpi = dpivalue)
plt.xscale('log')
plt.xlabel('$E_R$ (keV)')
plt.ylabel('PKA spectrum (Weighted) (a.u)')
plt.step(recEnbin, countpkaNRT, where = 'post')
plt.title(title)
plt.legend(elem)
plt.show()

plt.figure(4, dpi = dpivalue)
plt.xscale('log')
plt.xlabel('$E_R$ (keV)')
plt.ylabel('Cumulative Fraction of PKA less than $E_R$ (Weighted)')
for k in range(Ntargelm):
	ydata = []
	for j in range(n_recEn):
		if j > 0:
			ydata.append(pka_LTErNRT[j][k])
	plt.plot(xdata, ydata)	
plt.title(title)
plt.legend(elem)
plt.show()

ofile.close()

