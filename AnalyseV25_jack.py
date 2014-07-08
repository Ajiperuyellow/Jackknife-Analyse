#! /usr/bin/env python3.2
from numpy import *
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import math
from pylab import *
from scipy import *
from scipy.stats import *
import scipy as sy
from os import *
#from scipy import optimize
from scipy.optimize import curve_fit
#from numpy import polyfit
import numpy as N
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.figure
import codecs

#GIT Tracked
#CHECKIT

# This program runs through standard BAMPS-Hydro outputs and correlates time-dependent data
print("Go")
#Settings:

#1) Select the Column, in which you want to take the data out:
Timecolumn = 0
DatacolumnD1 = 15
DatacolumnD2 = 16
DatacolumnD3 = 17

#2) Select if you want to print out histograms and the data itself for each run:
histograms_on = False
currentplots_on = False
jackknife_output_on = True

#3) Data about the Runs:
run_group_name = 'RCONhX03NEW'
Number_of_runs = 20
timesteps = 10001
size_of_timestep = 0.001
Dimensions = 3

#4) Physical Data:
temperature_begin = 0.4
temperature_end = 0.4
temperature_inc = 0.1

efield = 0.0
efield_end = 0.0
efield_inc = 0.005

#5) Which Species are active and which are not:

Ups_Aktiv = True
Antiups_Aktiv = True
Downs_Aktiv = True
Antidowns_Aktiv = True
Strange_Aktiv = True
Antistrange_Aktiv = True
Gluons_Aktiv = True

#5) Name for the output-file:
beliebig = 'jackknifedev'

#6) Some Parameters:
ladung = [0.0, 2.0 / 3.0, -2.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0]
TT = [20000.0,2253.0,668.0,454849.0,232883.0,134771.0]#[34982.0,4373.0,1296.0,740.0] #    0.1 0.2 0.3 0.4  

#7) Settings about the Plots:

#Number of timesteps included in the Analysis!
kurzintervall=int(5000)
#number of X-labels. CAREFUL! kurzintervall/Anzahlxlabels HAS TO BE INTEGER!
Anzahlxlabels = 10
#Cutoff in the Beginning
Anfangscutoff = 0

#Printe:

print('ONLY EXPONENTIAL FIT')
print('*** Time-correlation/Current calculator V16 turbo ***')
print('Includes Testparticles!')
print('')

#Plotlabels
Positionxlabels = arange(0,kurzintervall,kurzintervall/Anzahlxlabels)
Arrayxlabels = arange(0,kurzintervall*size_of_timestep,kurzintervall*size_of_timestep/Anzahlxlabels)
kurzezeit=N.linspace(0,kurzintervall,kurzintervall)
kurzezeit_echtezeit=N.arange(0,kurzintervall*size_of_timestep,size_of_timestep)

#print('FILENAME: ' + filename_basis)
print('RUNS: ' + str(Number_of_runs))
print('timesteps: ' + str(timesteps))
print('Dimensions: ' + str(Dimensions))
print('')

#Diverses
Original_number_of_runs = Number_of_runs
Number_of_runs_times_dimensions = Number_of_runs*Dimensions
timesteps -= Anfangscutoff
GoForAll=True
filename_basis = ''
variancestring = ''
mittelwert = N.zeros(Number_of_runs_times_dimensions)
maxrunnumber = Number_of_runs + 1
Zeit_max = timesteps
zeit = range(0, Zeit_max)


#Begin with the Loops
b = temperature_begin
while b < (temperature_end + temperature_inc):
	a = efield
	while a < (efield_end + efield_inc):
		
		#Set Filenames
		filename_basis = run_group_name + "T" + str(b) + 'E' + str(a) + 'R'
		filename_basis = filename_basis.replace('.', '_')
		#Print Settings for this Loopstep
		print("*****************************************")
		print("file " + filename_basis)
		print("Efield = " + str(a))
		print("Temperatur = " + str(b))
		print("*****************************************")
		
		#Make/Change Folder
		mkdir('ANALYSIS_'+filename_basis+'_'+beliebig)
		chdir('ANALYSIS_'+filename_basis+'_'+beliebig)
		mkdir('Time_Histogram_'+filename_basis+'_'+beliebig)
		mkdir('Currents_'+filename_basis+'_'+beliebig)
		mkdir('Correlators_'+filename_basis+'_'+beliebig)
		chdir('../')	
		
		#Reset Variables
		variance_average = 0.0
		variance_std = 0.0
		Number_of_runs = Original_number_of_runs
		Number_of_runs_times_dimensions = Number_of_runs*Dimensions
		jackknife_string = ''

		#Initialize the arrays
		sum_of_corrs=N.zeros(Zeit_max)
		variance_array = N.zeros(Number_of_runs_times_dimensions)
		mean_corr_array = N.zeros(Zeit_max)
		std_corr_array = N.zeros(Zeit_max)
		corrfunction_array = N.zeros((Zeit_max, Number_of_runs_times_dimensions)) 
		corrfunction_array_manually = N.zeros((Zeit_max, Number_of_runs_times_dimensions)) 
		corrfunction_array_MASK = N.zeros((Zeit_max, Number_of_runs_times_dimensions))
		Array_aus_jackknife_steigungen = N.zeros(Number_of_runs_times_dimensions)
		
		#Loop through all runs
		runnumber = 1
		File_loaded = False
		No_file_count = 0

		while runnumber < maxrunnumber:
			#Print 
			if runnumber==50:
				print('Run: ' + str(runnumber))
			if runnumber==100:
				print('Run: ' + str(runnumber))
			if runnumber==150:
				print('Run: ' + str(runnumber))
			if runnumber==200:
				print('Run: ' + str(runnumber))
	
			#Initialize the arrays for THIS run
			vollerCorrelator = N.zeros(Zeit_max)
			stromcompleteD1 = N.zeros(Zeit_max)
			stromcompleteD2 = N.zeros(Zeit_max)
			stromcompleteD3 = N.zeros(Zeit_max)
			halberCorrelator = N.zeros(Zeit_max)
			richtigerCorrelator= N.zeros(Zeit_max)

			#####################################################
			#print('Run: ' + str(runnumber))
		
			#Read the files, numbered along this scheme name1.f11
			#File is in array of strings

			filename = filename_basis + str(runnumber) + '.f11'
			zeilen = []
			
			#Try-Clause because some files do not exist
			try:
				z = codecs.open(filename, mode='r')
				File_loaded = True
			except:
				print('File does not exist. Skip ' + filename)
				File_loaded = False
				No_file_count += 1
				#Set a Mask to this value that it is not subject to means, and so on...For all 3 dimensions
				corrfunction_array_MASK[:,runnumber-1]=1
				corrfunction_array_MASK[:,(runnumber+Number_of_runs)-1]=1
				corrfunction_array_MASK[:,(runnumber+(2*Number_of_runs))-1]=1
				
				
			#If the file Loaded, start the Analysis
			if File_loaded == True:
				for line in z.readlines():
					zeilen.append(line)
				z.close()
				
				#Set the Particle-Aktiv-Variable to Zero
				Aktiv = False
				
				#Loop through all lines for space-averaging of the current
				for line in zeilen:
					#Loops from top to bottom through lines and searches for keywords
					if not line.find('Particle-Type:1')==-1:
						Aktiv = Ups_Aktiv
						ch = 1
					if not line.find('Particle-Type:2')==-1:
						Aktiv = Antiups_Aktiv
						ch = 2
					if not line.find('Particle-Type:3')==-1:
						Aktiv = Downs_Aktiv
						ch = 3
					if not line.find('Particle-Type:4')==-1:
						Aktiv = Antidowns_Aktiv
						ch = 4
					if not line.find('Particle-Type:5')==-1:
						Aktiv = Strange_Aktiv
						ch = 5
					if not line.find('Particle-Type:6')==-1:
						Aktiv = Antistrange_Aktiv
						ch = 6						
						
					if Aktiv == True:
						#Split along Tabs
						col = line.split('\t')
						#Try-Clause because of Some lines wo data
						try:
							#Time sits in the first column
							zeitcounter = int(col[Timecolumn])
							zeitcounter -= Anfangscutoff
							#Data adds up in the Array
							if zeitcounter > -1:
								stromcompleteD1[zeitcounter] += float(col[DatacolumnD1])*ladung[ch]
								stromcompleteD2[zeitcounter] += float(col[DatacolumnD2])*ladung[ch]
								stromcompleteD3[zeitcounter] += float(col[DatacolumnD3])*ladung[ch]
							#print(str(runnumber) + '\t' + str(zeitcounter) + '\t' + str(stromcomplete[zeitcounter]))
						except:
							pass
				
				#The Average of D1 is saved in an array	
				mittelwert[runnumber-1]=N.mean(stromcompleteD1)
				
				#Plot now for THIS RUN Histograms or currents, D1
				#Plot Histograms
				if histograms_on == True:
					chdir('ANALYSIS_'+filename_basis+'_'+beliebig)
					chdir('Time_Histogram_'+filename_basis+'_'+beliebig)
					print('*** Plot Nx-Histogramm for this run ***')
					plt.hist(stromcomplete,100)
					plt.savefig(beliebig +'_HISTO-Current_' + '_'  + str(runnumber) + '.png', dpi=300, figsize=(8, 6))
					plt.clf()
					chdir('../../')
					
				#Plot Currents, D1
				if currentplots_on == True:
					correlationstring = ''
					for t in arange(Zeit_max):
						correlationstring += str(t) + '\t' + str(stromcompleteD1[t]) + '\n'
						t = t + 1
																	
					chdir('ANALYSIS_'+filename_basis+'_'+beliebig)
					chdir('Currents_'+filename_basis+'_'+beliebig)
					print('*** Plot Nx-Evolution for this run ***')
					#plt.plot(zeit, stromcomplete, linestyle='-', color='r', label=('N1'))
					#plt.xlabel('Time[steps]')
					#plt.ylabel('Current, MEAN: '+str(mittelwert[runnumber-1]))
					#plt.savefig('Part_'+str(typcount)+'_'+filename_basis+'_'+str(Number_of_runs)+'_runs_'+str(Zeit_max)+'_tsteps_'+beliebig+'.eps')
					#plt.savefig('Current' + '_T_'+ str(b) +'_'+ str(runnumber) + '.png', dpi=300, figsize=(8, 6))
					#plt.clf()
					
					#Write the output-file
					#chdir('ANALYSIS_'+filename_basis+'_'+beliebig)
					#chdir('Correlators_'+filename_basis+'_'+beliebig)
					f = codecs.open('Corr' + '_' + filename_basis + '_' + str(runnumber) + '_' + str(Zeit_max) + '_tsteps_' , 'w')
					f.write(correlationstring)
					f.close()
					chdir('../../')
					
				#CORRELATOR:
				#
				#1) Automatic correlator for all Dimensions
				#
				#D1
				vollerCorrelator = N.correlate(stromcompleteD1, stromcompleteD1, "full")
				halberCorrelator = N.array_split(vollerCorrelator, 2)
				richtigerCorrelator = N.flipud(halberCorrelator[0])
				N.concatenate((richtigerCorrelator.T, [0]))
				for t in arange(Zeit_max):
					richtigerCorrelator[t] /= float(Zeit_max - t)
					sum_of_corrs[t] += richtigerCorrelator[t]
					corrfunction_array[t,runnumber-1] =  richtigerCorrelator[t]
				variance_array[runnumber-1] = richtigerCorrelator[0]
				
				#D2
				vollerCorrelator = N.correlate(stromcompleteD2, stromcompleteD2, "full")
				halberCorrelator = N.array_split(vollerCorrelator, 2)
				richtigerCorrelator = N.flipud(halberCorrelator[0])
				N.concatenate((richtigerCorrelator.T, [0]))
				for t in arange(Zeit_max):
					richtigerCorrelator[t] /= float(Zeit_max - t)
					sum_of_corrs[t] += richtigerCorrelator[t]
					corrfunction_array[t,(runnumber+Number_of_runs)-1] =  richtigerCorrelator[t]
				variance_array[(runnumber+Number_of_runs)-1] = richtigerCorrelator[0]
				
				
				#D3
				vollerCorrelator = N.correlate(stromcompleteD3, stromcompleteD3, "full")
				halberCorrelator = N.array_split(vollerCorrelator, 2)
				richtigerCorrelator = N.flipud(halberCorrelator[0])
				N.concatenate((richtigerCorrelator.T, [0]))
				for t in arange(Zeit_max):
					richtigerCorrelator[t] /= float(Zeit_max - t)
					sum_of_corrs[t] += richtigerCorrelator[t]
					corrfunction_array[t,(runnumber+(2*Number_of_runs))-1] =  richtigerCorrelator[t]
				variance_array[(runnumber+(2*Number_of_runs))-1] = richtigerCorrelator[0]	
				
				#2) Manual Correlator
				#
				#Calculate the Variance manually if you want
				#variance = sum(pow(stromcomplete[i],2.0) for i in xrange(0, Zeit_max))/Zeit_max
				#
				#print('Calculate the Correlator manually')
				#print('---------------------------')		
				#t=0
				#while t < Correlator_maxtime:	
				#	corrfunction_array_manually[t,runnumber-1] =  sum(stromcomplete[i]*stromcomplete[i+t] for i in arange(Zeit_max-t))/(Zeit_max-t)
				#	t = t + 1
				
				
							
			runnumber = runnumber + 1
		###############################################################################
		###############################################################################	
		
		#JACKKNIFE
		
		#1) do the Analysis for all runs
		
		#Generate a clean array of corrfunctions, where  possible missing runs are masked away
		corrfunction_array_masked_clean = ma.masked_array(corrfunction_array, mask=corrfunction_array_MASK)*TT[int(b*10.0-1.0)]
		#Average the Correlator over all runs
		mean_corr_array = N.mean(corrfunction_array_masked_clean,axis=1)
		std_corr_array = N.std(corrfunction_array_masked_clean,axis=1)
		#Fit the Correlator to find the Exponential Decay
		def fitFunc(zeit, a, b):
			return a*N.exp(-1.0*zeit/b)
		#Cut the data off at kurzintervall
		daten_zum_fitten=mean_corr_array[0:kurzintervall]
		#Cut the errordata off at kurzintervall, divide by sqrt(N). N is now the actually used number of runs
		fehlerdaten_zum_fitten=std_corr_array[0:kurzintervall]/sqrt((Number_of_runs - No_file_count)*Dimensions)
		#DO THE LEAST_SQUARE FIT
		fitParams, fitCovariances = curve_fit(fitFunc, kurzezeit, daten_zum_fitten, p0=[1.0,1.0],sigma=fehlerdaten_zum_fitten)
		sigma = [sqrt(fitCovariances[0,0]), sqrt(fitCovariances[1,1])]
		#Save the value for the relaxation time for the full sample
		Full_sample_fit_value = fitParams[1]*size_of_timestep

		#2) do the Analysis for all runs but one, "Jackknifing""
		
		Array_aus_jackknife_steigungen = N.zeros(0)
		
		for jackknife_skip_number in arange(Number_of_runs_times_dimensions):
			#Number of runs used in each jackknife-Step
			Jackknife_number_of_runs = (Number_of_runs - No_file_count)*Dimensions - 1
			print("Jackknife_number_of_runs: " + str(Jackknife_number_of_runs))
			
			#Some runs do not exist. Leave them out. They are masked from the loop before.
			lasse_aus_weil_dieser_run_fehlt = False

			#So check if there is a mask on this special runnumber
			if corrfunction_array_MASK[:,jackknife_skip_number-1].all()==1:
				lasse_aus_weil_dieser_run_fehlt = True

			# if everything is allright, so the run exists, start with jackknifing
			if lasse_aus_weil_dieser_run_fehlt == False:
				
				#Mask this run. It is jackknifed away!
				corrfunction_array_MASK[:,jackknife_skip_number-1]=1

				#Generate a clean array of corrfunctions, where the jackknifing was done, and possible missing runs are also masked away
				corrfunction_array_masked_clean = ma.masked_array(corrfunction_array, mask=corrfunction_array_MASK)*TT[int(b*10.0-1.0)]

				# REDO the masking
				corrfunction_array_MASK[:,jackknife_skip_number-1]=0

				#Average the Correlator over all runs
				mean_corr_array = N.mean(corrfunction_array_masked_clean,axis=1)
				std_corr_array = N.std(corrfunction_array_masked_clean,axis=1)

				#Fit the Correlator to find the Exponential Decay
				def fitFunc(zeit, a, b):
					return a*N.exp(-1.0*zeit/b)

				#Cut the data off at kurzintervall
				daten_zum_fitten=mean_corr_array[0:kurzintervall]

				#Cut the errordata off at kurzintervall, divide by sqrt(N). N is now the actually used number of runs
				fehlerdaten_zum_fitten=std_corr_array[0:kurzintervall]/sqrt(Jackknife_number_of_runs)

				#DO THE LEAST_SQUARE FIT
				fitParams, fitCovariances = curve_fit(fitFunc, kurzezeit, daten_zum_fitten,p0=[1.0,1.0],sigma=fehlerdaten_zum_fitten)
				sigma = [sqrt(fitCovariances[0,0]), sqrt(fitCovariances[1,1])]

				#Save in Array
				N.append(Array_aus_jackknife_steigungen,fitParams[1]*size_of_timestep)

		print("Array_aus_jackknife_steigungen:")
		print Array_aus_jackknife_steigungen

		print("Summe:")
		print(sum((Array_aus_jackknife_steigungen-Full_sample_fit_value)**2)) 
		
		Jackknife_error = sqrt( (Jackknife_number_of_runs-1)/Jackknife_number_of_runs*sum((Array_aus_jackknife_steigungen-Full_sample_fit_value)**2) )
		print(" Jackknife-error estimation: " + str(Jackknife_error))


				
		#Update the number of runs
		Number_of_runs = Number_of_runs - No_file_count
		Number_of_runs_times_dimensions = Number_of_runs * Dimensions

		
		#Mask the Arrays to shield of Zero-Values because of missing runs
		corrfunction_array_masked_clean = ma.masked_array(corrfunction_array, mask=corrfunction_array_MASK)*TT[int(b*10.0-1.0)]


		#JACKKNIFE OUTPUT for Marc wagner's Code GEP'
		if jackknife_output_on:
			for s in arange(Number_of_runs_times_dimensions):
				for t in arange(Zeit_max):
					jackknife_string += str(t) + ' 0 0 ' + str(s) + ' ' + str(corrfunction_array_masked_clean[t,s]) + ' 0\n'
			print('\n jackknife string produced.')
			#Write the output-file
			chdir('ANALYSIS_'+filename_basis+'_'+beliebig)
			chdir('Correlators_'+filename_basis+'_'+beliebig)
			f = codecs.open(beliebig +'_JACKKNIFE_' + filename_basis + '_' + str(Number_of_runs_times_dimensions) + '_runs_' + str(Zeit_max) + '_tsteps_' , 'w')
			f.write(jackknife_string)
			f.close()
			chdir('../../')

		#Calculate the JACKKNIFE ERROR
		#print(Array_aus_jackknife_steigungen)
		#jackknife_error = sqrt((Number_of_runs_times_dimensions-1)*sum(pow((Array_aus_jackknife_steigungen[i]-fitParams[1]),2.0)/Number_of_runs_times_dimensions for i in arange(Number_of_runs_times_dimensions-1)))
		#print("old JACKKNIFE : " + str(jackknife_error))

		
		#Print out Information
		print("\n\n\n********* FINAL RESULTS *********\n\n\n")
		print(str(No_file_count) + ' files did NOT exist!!!' + '\n')
		
		print("1) VARIANCES")
		#REDUNDANT
		#print("This is redundant:!?")
		#print("==> Auto    Variance: 	" + str(sum_of_corrs[0]*TT[int(b*10.0-1.0)]/Number_of_runs_times_dimensions))
		
		#REDUNDANT
		#print("==> Average Variance: 	" + str(variance_average))
		#print("==> Uncertainty:      	" + str(variance_std/sqrt(Number_of_runs_times_dimensions)))
		
		print("==> Array Variance: 	" + str(mean_corr_array[0]))
		print("==> Array Uncertainty:	" + str(std_corr_array[0]/sqrt(Number_of_runs_times_dimensions)))

		print("==> Mean Current: 	" + str(N.mean(mittelwert)))
		print("==> Current Uncertainty:	" + str(N.std(mittelwert)/sqrt(Number_of_runs_times_dimensions)))		
		
	
		#Generate the output-Sring for the final information
		
		variancestring = '\n\n'
		variancestring += 'Variance + error + relaxation time + error:'
		variancestring += str(mean_corr_array[0]) + '\t' + str(std_corr_array[0]/sqrt(Number_of_runs_times_dimensions)) + '\t' + str(fitParams[1]*size_of_timestep) + '\t' + str(sigma[1]*size_of_timestep)+ '\n\n'
		#variancestring += str(b) + '\t' + str(variance_average*TT[int(b*10.0-1.0)]) + '\t' + str(variance_std*TT[int(b*10.0-1.0)]/sqrt(Number_of_runs_times_dimensions)) + '\n' + 'Mean Current =' + '\t' + str(N.mean(mittelwert)) + '\t' + '+-' + '\t' + str(N.std(mittelwert)/sqrt(Number_of_runs)) + '\n'
		variancestring += '*** STATISTICAL TEST on Gaussian property of the Variance *** \n'+ '\n'
		variancestring +='1) ANDERSON-Test: \n'
		variancestring +='A2 = ' + str(A2) + '\n'
		variancestring +=str(sig) + '\n'
		variancestring +=str(crit) + '\n'+ '\n'
		variancestring +='2) KS-Test: \n'
		variancestring +='p-Value = ' + str(pD) + '\n'
		variancestring +='W-Value = ' + str(D) + '\n'+ '\n'
		variancestring +='3) Shapiro-Test:\n'
		variancestring +='p-Value = ' + str(pW) + '\n'
		variancestring +='W-Value = ' + str(W) + '\n'	+ '\n'
		variancestring += '*** STATISTICAL TEST on Gaussian property of the Mean Current *** \n'+ '\n'
		variancestring +='1) ANDERSON-Test: \n'
		variancestring +='A2 = ' + str(StromA2) + '\n'
		variancestring +=str(Stromsig) + '\n'
		variancestring +=str(Stromcrit) + '\n'+ '\n'
		variancestring +='2) KS-Test: \n'
		variancestring +='p-Value = ' + str(StrompD) + '\n'
		variancestring +='W-Value = ' + str(StromD) + '\n'+ '\n'
		variancestring +='3) Shapiro-Test:\n'
		variancestring +='p-Value = ' + str(StrompW) + '\n'
		variancestring +='W-Value = ' + str(StromW) + '\n'	
		variancestring +=  '\n\n' +  '*** RELAXATION-TIME: ***\n' +  str(fitParams[1]*size_of_timestep) + ' +- ' + str(sigma[1]*size_of_timestep)		
		
		#Generate the output-String for the Full Correlator
		
		superstring = ''
		for t in arange(Zeit_max):
			#sum_of_corrs[t] /= Number_of_runs
			superstring += str(t) + '\t' + str(mean_corr_array[t])+ '\t' + str(std_corr_array[t]/sqrt(Number_of_runs_times_dimensions))  + '\n'
			t = t + 1

		#Print Information
		print('Do the All-Particle Plots')
		
		#Write the output-file
		chdir('ANALYSIS_'+filename_basis+'_'+beliebig)
		chdir('Correlators_'+filename_basis+'_'+beliebig)
		f = codecs.open(beliebig +'Corr_' + '_' + filename_basis + '_' + str(Number_of_runs_times_dimensions) + '_runs_' + str(Zeit_max) + '_tsteps_' , 'w')
		f.write(superstring)
		f.close()

		#Plot and save the figure
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.xaxis.grid
		ax.yaxis
		figure.autolayout = True
		
		Atex=r'$'+str(round(fitParams[0],4))+' \pm ' + str(round(sigma[0],4)) + ' fm^{-6}$'
		tautex=r'$'+str(round(fitParams[1]*size_of_timestep,5))+' \pm ' + str(round(sigma[1]*size_of_timestep,4)) + ' fm$'
		#ctex = r'$ c = '+ str(round(fitParams[2],4)) + ' \pm ' + str(round(sigma[2],4)) + ' fm^{-6}$'
		
		tex = r'$C(0)= $' + Atex + '\n ' + r'$\tau= $' + tautex + '\n' #+ ctex
		#tex = r'$C(0)= $' + Atex + '\n ' + r'$\tau= $' + tautex + '\n' + ctex
		#ax.text(fitParams[1]/2.0, fitParams[0]+fitParams[2], tex, fontsize=15, va='bottom')
		ax.text(fitParams[1]/2.0, fitParams[0], tex, fontsize=15, va='bottom')
		
		ax.set_xlabel(r'Time $[fm]$')
		ax.tick_params(axis='x', pad=-0.7)
		
		#ax.errorbar(x, y, yerr=[yerr_lower, 2*yerr], xerr=xerr,fmt='o', ecolor='g', capthick=2)
		plt.grid(b=True, which='major')
		#plt.plot(kurzezeit, fitFunc(kurzezeit, fitParams[0], fitParams[1], fitParams[2]))
		plt.plot(kurzezeit, fitFunc(kurzezeit, fitParams[0], fitParams[1]))#,zeit, fitFunc(t, fitParams[0] + sigma[0], fitParams[1] - sigma[1], fitParams[2] + sigma[2]),zeit, fitFunc(t, fitParams[0] - sigma[0], fitParams[1] + sigma[1], fitParams[2] - sigma[2]))
		plt.errorbar(kurzezeit, daten_zum_fitten, linestyle='-', color='r', label=('Corr'), yerr=fehlerdaten_zum_fitten, ecolor='g')

		plt.ylabel(r'$\sigma/T$')
		plt.title(r'Time-Correlator of $N^1$')
		#plt.legend(loc=2,shadow=True)
		#plt.yscale('log')
		#plt.xscale('log')
		#plt.ylim(0.0,5)
		plt.xlim(0,kurzintervall)
		#xticks = Npy.arange(0.6,2.57,0.18)
		#yticks = Npy.arange(0.007,0.1,0.01)
		#plt.xaxis.set_ticks( xticks )
		#plt.yaxis.set_ticks( yticks )
		#plt.yticks(yticks,yticks)
		ax.set_xticks(Positionxlabels)
		ax.set_xticklabels(Arrayxlabels,rotation=45)
		
		plt.ylabel(r'Correlator C(t) $[fm^{-6}]$')
		ax.bottom = 0.25
		
		#plt.savefig('Part_'+str(typcount)+'_'+filename_basis+'_'+str(Number_of_runs)+'_runs_'+str(Zeit_max)+'_tsteps_'+beliebig+'.eps')
		plt.savefig(beliebig +'Corr_' + '_' + filename_basis + '_' + str(Number_of_runs_times_dimensions) + '_runs_' + str(Zeit_max) + '_tsteps_' +  '.png', dpi=300, figsize=(8, 6))
		plt.clf()
		f = codecs.open( beliebig+'VARIANCE' + '_' + filename_basis + '_' + str(Number_of_runs_times_dimensions) + '_runs_' + str(Zeit_max) + '_tsteps_' , 'w')
		f.write(variancestring)
		f.close()
		#Change the directory back
		chdir('../../')
		#Increase the Variables
		a = a + efield_inc
	b += temperature_inc





    ##########################################################################################################################################################################
    ##########################################################################################################################################################################
    ##########################################################################################################################################################################



print('Finish')

