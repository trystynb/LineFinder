#!/usr/bin/python2.7

""" 
TB_LINEFINDER.PY
v1.2.1 (5/02/2016)
	Previous version was not loading limits properly from log. This is now fixed.
v1.2 (11/01/2016)
	MAJOR revisions are:
		-Changed element selection in LINEADDER from radiobuttons to checkbox
			(fixes unselecting elements)
		-Changed log-writing/editing to a GUI window rather than through terminal window (UPDATELOG class)
			-Can now left/right clcik on spectrum display to set velocity limits!
		-Fixed behaviour of log-editing drop down menus (they now appear to update when a line/system/ion is removed)
			-NOTE: This doesn't fix the bug if there isn't a redshift/line/ion in the log file, but
				it at least will notify the user that the selected redshift is "<Z>", etc
		-Tutorial mode is now off by default!
	MINOR revisions not noted, but minor code adjustments to plots, code, etc.
	

V1.1 (10/07/2015)
	-Added a file browser button for selecting files.
	-Added a tutorial mode that can be turned on/off in the main menu.

Author: Trystyn Berg (trystynb@uvic.ca)

Stand-alone python application to identify spectral lines and save selected line list.


PYTHON PACKAGE REQUIREMENTS:
	TB_LINEFINDER was developed with:
		-Python Version 2.7
		-Matplotlib Version 1.4.1
			-Need matplotlib.backends.backend_tkagg
		-Tkinter Revision: 81008
		-Numpy
	Other verisons of above packages may work...

INPUT FILE REQUIREMENTS
	-ASCII file containing spectrum for plotting.
		Each line in file must be in the format (whitespace delimeted):
			WAVELENGTH FLUX
		Any additional columns after first two are ignored.
		MUST CONTAIN THE SAME NUMBER OF COLUMNS IN EACH LINE (unless commented)
		Can use '#' at the start of any line to remove comment

	-ASCII file containing a line list for identification.
		Each line in the file must contain information on each spectral line
		of interest in the format (whitespace delimited):
			REST_WAVELENGTH ION SHORT_WL OSCILLATOR_STRENGTH

			REST_WAVELENGTH - The rest wavelength of the spectral line (e.g. 1215.6701)
			ION - The species (including ionization state) of the spectral line (e.g. HI)
			SHORT_WL - A short identifier of the wavelength (e.g. 1215)
				Needed to descrimante between absorption lines of a given ION
			OSCILLATOR_STRENGTH - The oscillator strength of the feature (e.g. 0.4164)
				This is only used for display purposes.

			EXAMPLE:
			1215.6701 HI 1215 0.4164

		If Ly-alpha not present, TB_LINEFINDER will add it by default
			(Morton 2003 values, as in example above)

		MUST CONTAIN THE SAME NUMBER OF COLUMNS IN EACH LINE (unless commented)

		Can use '#' at the start of any line to remove or comment

OPTIONAL INPUT
	-TB_LINEFINDER log file (see OUTPUT FILE below for description)
OUTPUT FILE
	-TB_LINEFINDER will generate a log file with the line list.
		The output is a semicolon (;) delimited
		Each line represents an spectral line identified in the format:
			Z; ION; LINE; FLAG; VMIN; VMAX; NOTES; COLOUR;
			Z - Redshift (to 5 decimal places) of line identified
			ION - The Species of the line identified (from input linelist)
			LINE - The SHORT_WL identifier from the input linelist
			FLAG - A integer associated with the quality of the identified line
			VMIN - The bluemost velocity (km/s) of the line profile
			VMAX - The redmost velocity (km/s) of the line profile
			NOTES - A string with any notes the user inputs
			COLOUR - The matplotlib colour to denote feature in TB_LINEFINDER
			
		The FLAG notation is a binary format with the following option
			0 - No good/skip
			1 - OK
			2 - Blend
			4 - Upper Limit/Non-detection
			8 - Lower Limit/Saturated


		When the output logfile is being read, any line beginning with '#'
		will be ignored

		For the purpose of future reference, the first two lines of the
		logfile contain the input filenames in the format:

		#!Line List: <INPUT LINELIST ASCII FILENAME>
		#!Spectrum: <INPUT SPECTRUM ASCII FILENAME>



KNOWN ISSUES:
	Often MATPLOTLIB will change some features, which affect how the user can
	interact with the application. If problem, please check your matplotlib version
"""

##############################
###INTIALIZATION OF PROGRAM###
##############################

#Import necessary Matplotlib packages
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
#Import Tkinter for GUI interface
import Tkinter, tkFileDialog, tkMessageBox
#import other basic pacakages
import numpy as np
import math
import os,sys


#Default input/output files necessary for running
initinspec='test.ascii'
initllist='linelist.lst'
initlog='linefinder.log'


#Command line inputs (Can be changed or set within GUI)
inputs=sys.argv
if len(sys.argv)==4:
	initinspec=sys.argv[1]
	initllist=sys.argv[2]
	initlog=sys.argv[3]

#List of colours to use (must be Matplotlib colours)
mplcolours=['r','k','b','c','m','g']
#DEBUG set to TRUE enables print statements throughout code for debugging purposes
debug=False
#debug=True


#USETUTORIAL is a boolean to control message boxes that help use the software.
#Thic can be turned on/off using the menu in the main Linefinder window.
#usetutorial=True
usetutorial=False

def velspec(line, wvlngth, vmin,vmax, spectrum):
	""" 
	VELSPEC takes a spectrum and converts the wavelength scale
	into a velocity scale (in km/s) based on an input wavelength

	call - TMPVEL,TMPSPEC=velspec(LINE,WVLNGTH,VMIN,VMAX,SPECTRUM)

	INPUT VARIABLES:
		LINE - The wavelength to zero the velocity scale to (float)
		WVLGNTH - A NUMPY array with the wavelength information of the spectrum
		VMIN - The minimum velocity of the returned velocity scale (float)
		VMAX - The maximum velocity of the returned velocity scale (float)
		SPECTRUM - A NUMPY array with the flux information of the spectrum

	RETURNS:
		TMPVEL - A NUMPY array of the Velocity-convereted scale (between VMIN and VMAX)
		TMPSPEC - A slice of the input SPECTRUM, with each element corresponding to the
			flux at the same element in the outputted TMPVEL

	NOTES:
		If there is no coverage of the velocity specturm between VMIN and VMAX,
		TMPSPEC will contain an array of zeros. If this occurs, a warning will
		be printed to the screen.

	"""
	c=3.0E5# speed of light (km/s)
	vel=np.array((wvlngth-line)/line*c)#Convert to a velocity
	#Check to see that VEL is within the specified VMIN/VMAX  range
	if min(vel)<vmin and max(vel)>vmax:
		ind1=(min(np.where(vel>vmin)[0]))#Find the index closest to vmin
		ind2=(max(np.where(vel<vmax)[0]))#Find the index closest to vmax
		tmpvel=vel[ind1:ind2]#Generate the velocity (x) array of intere$
		tmpspec=spectrum[ind1:ind2]#Do the Same for the relative flux(y)
		return tmpvel, tmpspec#Return slice of spectrum
	#If not, return zero arrays and print warning
	else:
		print 'VELSPEC: Line not in velocity range'
		return vel, np.zeros(len(vel))
def specfits(infile):	
	""" 
	SPECFITS reads the input spectrum ASCII file and returns the data in two NUMPY arrays

	Call - WVLNGTH,SPECTRUM=specfits(INFILE)

	INPUT VARIABLE:
		INFILE - The input spectrum file.
			Each line in file must be in the format (whitespace delimeted):
				WAVELENGTH FLUX
			Any additional columns after first two are ignored.
			MUST CONTAIN THE SAME NUMBER OF COLUMNS IN EACH LINE (unless commented)
			Can use '#' at the start of any line to remove comment

	OUTPUT VARIABLES:
		WVLNGTH - A NUMPY array with the wavelength (first column) of the input spectrum
		SPECTRUM - A NUMPY array with the flux (second column) of the input spectrum


	NOTES:
		If no input spectrum file is found, WVLNGTH and SPECTRUM
		will be NUMPY arrays full of zeros. A warning will be printed
		to the screen.

	"""
	#Check to see if file exists
	if os.path.isfile(infile):
		#REad spectrum file. All values must be floats
		#Will ignore lines with '#' flag
		#First column must be wavelength
		#Second column must be flux
		data=np.genfromtxt(infile,dtype=type(0.00),comments='#')
		wvlngth=data[:,0]
		spectrum=data[:,1]
		#Return spectrum data
		return wvlngth,spectrum
	#If file doesn't exist
	else:
		#Print warning
		print "WARNING: No spectrum file found"
		#Return zero arrays
		return np.zeros(100),np.zeros(100)



def IsFloat(str):
	""" 
	ISFLOAT checks if a given string is a float.
	
	Call - OUTPUT=IsFloat(STR)
	
	INPUT VARIABLE:
		STR - A string for checking whether a float or not
	OUTPUT VARIABLE:
		OUTPUT - A boolean value corresponding to whether STR
		is a float (returns TRUE) or not (returns FALSE)
	
	"""
	#See if STR is a float
	try:
		float(str)
	#If test fails, return FALSE
	except Exception:
		return False
	#If passed the test, return TRUE
	return True


def WriteLog(log,logfile,llistfile,fits):
	""" 
	WRITELOG writes the TB_LINEFINDER output logfile.
	
	call WriteLog(LOG,LOGFILE,LLISTFILE,FITS)
	
	INPUT VARIABLES:
		LOG - A dictionary with the log information (see LOADLOG for structure)
		LOGFILE - A string with the desired filename for the output logfile
		LLISTFILE - A string with the input linelist filename
		FITS - A string with the input spectrum filename

	NOTES:
		If no logfile is provided, it will not save log and display warning
	
	"""
	#Check to see if LOGFILE string has at least one character
	if len(logfile)>1:
		#Open LOGFILE buffer for writing
		f=open(logfile,'w')
		#Write the input files used to generate TB_LINEFINDER logfile
		f.write('#!Line List: %s\n'%llistfile)
		f.write('#!Spectrum: %s\n'%fits)
		#Write column labels for user reference
		f.write('#z;\t\tIon;\tline;\tflag;\tvmin;\tvmax;\tNotes;\tcolour\n')

		#Loop through the LOG dictionary (redshift, ion/species, line identifier)
		for z in log['zs']:
			for ion in log['ions']:
				for line in log['lines'][ion]:
					#check if entry (keyed by z,ion,line) exists in log
					if (z,ion,line) in log:
						#If so, write to log file in format:
						#Z; ION; LINE; FLAG; VMIN; VMAX; NOTES; COLOUR;
						f.write('%s;\t'%z)
						f.write('%s;\t'%ion)
						f.write('%s;\t'%line)
						f.write('%s;\t'%log[z,ion,line,'flag'])
						f.write('%s;\t%s;\t'%log[z,ion,line,'vel'])
						f.write('%s;\t'%log[z,ion,line])
						f.write('%s;\n'%log[z,ion,line,'colour'])
		#close Logfile writing buffer
		f.close()
		#Print stuff informing user of WRITELOG's success (including contents
		#of file
		print "Wrote logfile: %s"%logfile
		print "Contents:\n"
		print "***********\n"
		os.system('cat %s'%logfile)
		print "\n***********"
	#If LOGFILE doesn't have a single character, cannot save. Print Warning
	else: print "No log file provided. Did not Save."
def LoadLog(logfile):
	""" 
	LOADLOG loads the log file.


	
	Call: LoadLog(LOGFILE)

	INPUTS:
		LOGFILE - String containing filename out TB_LINEFINDER output logfile.
		The format of the log file is as follows:
			For the purpose fo future reference, the first two lines of the
			logfile contain the input filenames in the format:
				#!Line List: <INPUT LINELIST ASCII FILENAME>
				#!Spectrum: <INPUT SPECTRUM ASCII FILENAME>

			The file is a semicolon (;) delimited
			Each line represents an spectral line identified in the format:
				Z; ION; LINE; FLAG; VMIN; VMAX; NOTES; COLOUR;
				Z - Redshift (to 5 decimal places) of line identified
				ION - The Species of the line identified (from input linelist)
				LINE - The SHORT_WL identifier from the input linelist
				FLAG - A integer associated with the quality of the identified line
				VMIN - The bluemost velocity (km/s) of the line profile
				VMAX - The redmost velocity (km/s) of the line profile
				NOTES - A string with any notes the user inputs
				COLOUR - The matplotlib colour to denote feature in TB_LINEFINDER
			
			The FLAG notation is a binary format with the following option
				0 - No good/skip
				1 - OK
				2 - Blend
				4 - Upper Limit/Non-detection
				8 - Lower Limit/Saturated

			When the output logfile is being read, any line beginning with '#'
			will be ignored

	OUTPUT:
		LOG - A dictionary with all the log information.
		The dictionary is keyed in the following way

		LOG['zs'] - A list of redshifts for each system

		LOG['ions'] - A list of all ion species (e.g. HI, CIV) within LOG.
			All ion species should be found in the linelist input to TB_LINEFINDER

		LOG['lines']- A dictionary of all the lines within LOG. This dictionary
			is keyed by all the ion species in LOG['ions'].
			For each ION in LOG['ions']:
				LOG['lines'][ION] - A list of all the lines (identified by the
				SHORT_WL identifier in the linelist input to TB_LINEFINDER

		LOG[Z,ION,LINE] - COntains the NOTES information in the LOG dictionary. It
			is keyed by (Z,ION,LINE), which are found in the lists of
			LOG['zs'], LOG['ions'], and LOG['lines'][ION] (respectively)

		LOG[Z,ION,LINE,'flag'] - COntains the FLAG information in the LOG dictionary. It
			is keyed by (Z,ION,LINE,'flag'), which are found in the lists of
			LOG['zs'], LOG['ions'], and LOG['lines'][ION] (respectively)

		LOG[Z,ION,LINE,'colour'] - COntains the COLOUR information in the LOG dictionary. It
			is keyed by (Z,ION,LINE,'colour'), which are found in the lists of
			LOG['zs'], LOG['ions'], and LOG['lines'][ION] (respectively)

		LOG[Z,ION,LINE,'vel'] - COntains the tuple (VMIN,VMAX) from the LOG dictionary. It
			is keyed by (Z,ION,LINE,'vel'), which are found in the lists of
			LOG['zs'], LOG['ions'], and LOG['lines'][ION] (respectively)

	NOTES:
		If no LOGFILE string is found in the current directory, a blank 
		LOG dictionary is loaded, with LOG['ion'], LOG['zs'], and LOG['lines']
		generated. The file will be created at saving the log.

	"""
	#THe main LOG dictionary
	log={}
	#Key LOG withe the redshift and ion lists
	log['zs']=[]
	log['ions']=[]
	#Key LOG with a sub-dictionary of the spectra lines
	log['lines']={}
	#Check if the LOGFILE all ready exists, and read it.
	if os.path.isfile(logfile):
		print "Loading logfile: %s"%logfile
		#Read in the semicolon delimted file
		data=np.genfromtxt(logfile,delimiter=';',dtype=type('str'),comments='#')
		if len(data)>0:
			#For each line in the file
			for ii in range(len(data)):
				#REad in each column, remove whitespace around the string,
				#and format the data.
				z=data[ii][0].strip()
				ion=data[ii][1].strip()
				line=data[ii][2].strip()
				#If Z isn't in 'zs' list, add it.
				if z not in log['zs']: log['zs'].append(z)
				#If ION isn't in 'ions' list, add it.
				#Also, need to add ION as a key to 'lines' dictionary
				if ion not in log['ions']:
					log['ions'].append(ion)
					log['lines'][ion]=[]
				#If LINE isn't in ['lines'][ION], include it
				if line not in log['lines'][ion]: log['lines'][ion].append(line)
				log[z,ion,line,'flag']=data[ii][3].strip()
				log[z,ion,line,'vel']=data[ii][4].strip(), data[ii][5].strip()
				log[z,ion,line]=data[ii][6].strip()#Add NOTE
				log[z,ion,line,'colour']=data[ii][7].strip()
	#If no LOGFILE is found, inform the user that is the case.
	else:
		print "Logfile not found at startup. Will create %s on exit"%logfile
	#return LOG
	return log


class UpdateLog:
	""" 
	CLASS UPDATELOG - The GUI environment for editing and updating the LOG file/dictionary for a given
		absorption feature. A plot of the absorption line will be displayed, along with a window to
		set the log for the absortption feature

		This will allow the user to set the velocity bounds of the absroption
		feature (left/right clicking on the display window), the flags of the line, any notes/comments,
		and the colour for the LINEFINDER display. A

	Call - UL=UpdateLog(LOG,Z,ION,LINE,WVL,IWVLNGTH,ISPECTRUM)

	INPUTS:
		LOG - The log dictionary 
		Z - The redshift (float) of the system of interest
		ION - The name of the ion in the log file that is to be edited
		LINE - The wavelength identifier of the line to be edited
		WVL - The real wavelength (float) of the absorption line (from LLIST?)
		IWVLNGTH - The array of wavelengths of the input spectrum
		ISPECTRUM - The array of fluxes of the input spectrum
	ATTRIBUTES:
		UL.LOG - The log dictionary
		UL.ZSTR - The String of redshift
		UL.ION - The ion of the absroption feature
		UL.LINE - The wavelength identifier of the absorption feature
		UL.ULMASTER - The TKinter window master
		UL.ULFRAME - The frame (embedded in UL.ULMASTER) to contain the velocity plot
		UL.OKFLAG - Boolean keeping track if the line is OK to use
		UL.RBOK - A Checkbox button that sets UL.OKFLAG
		UL.BLENDFLAG - Boolean keeping track if the line is blended
		UL.RBBLEND - A checkbox button that sets UL.BLENDFLAG
		UL.UPPERFLAG - Boolean keeping track if the line should be classified as upper limit
		UL.RBUPPER - A checkbox that sets UL.UPPERFLAG
		UL.LOWERFLAG - Boolean keeping track if the line is a lower limit
		UL.RBLOWER - A Checkbox that sets UL.LOWERFLAG
		UL.COLOUR - String of matplotlib colour to flag the feature with on main LINEFINDER display window
		UL.NOTES - The string of notes associated with the absorption feature
		UL.VMIN - The minimum (blueward) velocity tto start absorption feature
		UL.VMAX - The maximum (redward) velocity to end absorption feature
		UL.FIG - Matplotlib figure (embedded in UL.ULFRAME) to diosplay velocity profile of line
		UL.AX - Matplotlib axes instance of velocity profile
		UL.VMINLINE - The vertical line object to signify UL.VMIN on UL.AX plot
		UL.VMAXLINE - The vertical line object to signify UL.VMAX on UL.AX plot
		UL.CANVAS - The TKinter canvas in UL.ULFRAME for displaying plots.
		UL.DRAWVMIN(REDRAW=TRUE) - (Re)draws the UL.VMIN line (UL.VMINLINE) on the UL.AX plot.
			Set REDRAW=FALSE if it is the first instance.
		UL.DRAWVMAX(REDRAW=TRUE) - (Re)draws the UL.VMAX line (UL.VMAXLINE) on the UL.AX plot.
			Set REDRAW=FALSE if it is the first instance.			
		UL.ONSETVELS(EVENT) - Takes the mouse click event, and gets the x-coordinate (velocity).
			If the left/right mouse button, UL.ONSETVELS fills in the entry for UL.VMIN/UL.VMAX
			value, and redraws the UL.VMINLINE/UL.VMAXLINE (respectively) 
		UL.ONSAVEBUTTON - Takes all the information (flags, notes, colours, velocity limits) for the
			absrotpion profile, saves it to UL.LOG, and destroys the UPDATELOG widget.

	NOTES:
		-The left mouse button shoudl set UL.VMIN on the window
		-The right mouse button sets UL.VMAX.
		-The flag coding is:
			0 - No good/skip
			1 - OK (i.e. use the line)
			2 - Blend
			4 - Upper Limit/Non-detection
			8 - Lower Limit/Saturated

	"""
	def __init__(self,log, z,ion,line,wvl, iwvlngth, ispectrum):
		#Display the tutorial message for how to use UPDATELOG
		if usetutorial: tkMessageBox.showinfo("Help Message", "Input parameters associated with the absroption line."+ \
			"Left/right clicking will set velocity bounds of absroption feature.")
		#Set the input information as attributes
		self.log=log
		self.zstr='%.5f'%z
		self.ion=ion
		self.line=line
		#Define the TK window
		self.ULmaster=Tkinter.Toplevel()

		#Create a frame for the velocity plot of the absroption line
		self.ULFrame=Tkinter.Frame(self.ULmaster, width=80,height=40)
		self.ULFrame.grid(row=0,column=0,columnspan=4)
		#Set the name of the TK window
		self.ULmaster.wm_title("Edit Log for z=%s (%s,%s)"%(self.zstr,self.ion,self.line))

		#Initialize the variables in the log dictionary if they don't exist all ready
		if (self.zstr,self.ion,self.line,'flag') not in self.log: self.log[self.zstr,self.ion,self.line,'flag']='0'
		if (self.zstr,self.ion,self.line,'colour') not in self.log: self.log[self.zstr,self.ion,self.line,'colour']='k'
		if (self.zstr,self.ion,self.line) not in self.log: self.log[self.zstr,self.ion,self.line]=''
		if (self.zstr,self.ion,self.line,'vel') not in self.log: self.log[self.zstr,self.ion,self.line,'vel']='-50','50'


		#Initialize different flags as boolean values 
		#Get the flag information from the log
		flags=int(self.log[self.zstr,self.ion,self.line,'flag'])
		#If flag contains an 8, it is a lower limit
		initLower=(flags>=8)
		#remove the 8 from the flag if necessary
		if initLower: flags-=8
		#If flag contains 4, it is an upper limit (remove if necessary)
		initUpper=(flags>=4)
		if initUpper: flags-=4
		#If flag contains a 2, it is a blend (remove again...)
		initBlend=(flags>=2)
		if initBlend: flags-=2
		#If flag is 1, it is OK to use!
		initOK=(flags==1)


		#Set up radio buttons for flags:
		#Line is OK (i.e. good to use)
		self.OKFlag=Tkinter.BooleanVar()
		self.OKFlag.set(initOK)
		self.RBOK=Tkinter.Checkbutton(self.ULmaster,text='Use line',variable=self.OKFlag,state='active')
		self.RBOK.grid(column=0,row=1)

		#Line is blended
		self.BlendFlag=Tkinter.BooleanVar()
		self.BlendFlag.set(initBlend)
		self.RBBlend=Tkinter.Checkbutton(self.ULmaster,text='Blend',variable=self.BlendFlag,state='active')
		self.RBBlend.grid(column=1,row=1)

		#Line should be used as an upper limit
		self.UpperFlag=Tkinter.BooleanVar()
		self.UpperFlag.set(initUpper)
		self.RBUpper=Tkinter.Checkbutton(self.ULmaster,text='Upper limit/non-detection',variable=self.UpperFlag,state='active')
		self.RBUpper.grid(column=2,row=1)

		#Line should be used as a lower limit
		self.LowerFlag=Tkinter.BooleanVar()
		self.LowerFlag.set(initLower)
		self.RBLower=Tkinter.Checkbutton(self.ULmaster,text='Lower limit/saturated',variable=self.LowerFlag,state='active')
		self.RBLower.grid(column=3,row=1)


		#Set-up drop down menu for colour coding
                labelCol=Tkinter.StringVar()#Name of string in label
                labelcol=Tkinter.Label(self.ULmaster,textvariable=labelCol,\
                        anchor="w",fg="black")
                labelcol.grid(column=0, row=2, columnspan=2, sticky='EW')
                labelCol.set(u"Colour for display:")#Initial value of Variable
		#SELF.COLOUR will be the colour of the line
		self.Colour=Tkinter.StringVar()
		self.Colour.set(self.log[self.zstr,self.ion,self.line,'colour'])
		CMenu=apply(Tkinter.OptionMenu, (self.ULmaster,self.Colour,)+tuple(mplcolours))
		CMenu.grid(column=2,row=2,columnspan=2,sticky='EW')

		#Set-up field for entering comments on the line profile
                labelNote=Tkinter.StringVar()#Name of string in label
                labelnote=Tkinter.Label(self.ULmaster,textvariable=labelNote,\
                        anchor="w",fg="black")
                labelnote.grid(column=0, row=3, sticky='EW')
                labelNote.set(u"Comments:")#Initial value of Variable
		#SELF.NOTES will have the notes for the line
		self.Notes=Tkinter.StringVar()
		entrynote=Tkinter.Entry(self.ULmaster,textvariable=self.Notes)
		entrynote.grid(column=1,row=3,columnspan=3,sticky='EW')
		self.Notes.set(self.log[self.zstr,self.ion,self.line])


		#Set-up field for entering bounding velocities of absorption
                labelVmin=Tkinter.StringVar()#Name of string in label
                labelvmin=Tkinter.Label(self.ULmaster,textvariable=labelVmin,\
                        anchor="w",fg="black")
                labelvmin.grid(column=0, row=4,sticky='EW')
                labelVmin.set(u"Blue velocity bound:")#Initial value of Variable
		self.Vmin=Tkinter.DoubleVar()
		entryvmin=Tkinter.Entry(self.ULmaster,textvariable=self.Vmin)
		entryvmin.grid(column=1,row=4,sticky='EW')
		self.Vmin.set(float(self.log[self.zstr,self.ion,self.line,'vel'][0]))

		#Set-up field for entering bounding velocities of absorption
                labelVmax=Tkinter.StringVar()#Name of string in label
                labelvmax=Tkinter.Label(self.ULmaster,textvariable=labelVmax,\
                        anchor="w",fg="black")
                labelvmax.grid(column=2, row=4,sticky='EW')
                labelVmax.set(u"Red velocity bound:")#Initial value of Variable
		self.Vmax=Tkinter.DoubleVar()
		entryvmax=Tkinter.Entry(self.ULmaster,textvariable=self.Vmax)
		entryvmax.grid(column=3,row=4,sticky='EW')
		self.Vmax.set(float(self.log[self.zstr,self.ion,self.line,'vel'][1]))
		
		#Set up the save and exit button
		SaveQuitButton=Tkinter.Button(self.ULmaster,text="Save Log & Exit", command=self.onSaveButton)
		SaveQuitButton.grid(column=0,row=5,columnspan=4, sticky='EW')


		#Set up Plot window
		#Velocity limit fields
		vmin=-1000.0
		vmax=1000.0
		#WVL should be Redshifted 
		wvl=wvl*(1.0+z)
		#get velocity profile of the line from the spectrum
		tmpvel,tmpspec=velspec(wvl, iwvlngth, vmin,vmax, ispectrum)
		self.fig=plt.figure(figsize=(6,6))
		self.fig.subplots_adjust(hspace=2.0,wspace=2.0)
		self.ax=plt.subplot(1,1,1)
		self.ax.plot(tmpvel,tmpspec,'k',drawstyle='steps')
		self.ax.set_xlabel(r'Relative velocity (km s$^{-1}$)')
		self.ax.set_ylabel(r'Relative intensity')
		self.Canvas=FigureCanvasTkAgg(self.fig,self.ULFrame)
		self.Canvas._tkcanvas.grid(row=0,column=0)
		self.DrawVmin(redraw=False)
		self.DrawVmax(redraw=False)
		self.fig.canvas.mpl_connect("button_press_event",self.OnSetVels)
		self.Canvas.start_event_loop(0)

	#Functions to draw velocity limits on UL.AX instance. One is for each limit
	def DrawVmin(self,redraw=True):
		#REDRAW defines whether the line has to be redrawn. If not, it needs to create UL.VMINLINE 
		ymin,ymax=self.ax.get_ylim()
		if redraw: self.ax.lines.remove(self.vminline[0])
		self.vminline=self.ax.plot([self.Vmin.get(),self.Vmin.get()],[ymin,ymax],'--r')
		self.Canvas.draw()
	def DrawVmax(self,redraw=True):
		#REDRAW defines whether the line has to be redrawn. If not, it needs to create UL.VMAXLINE 
		ymin,ymax=self.ax.get_ylim()
		if redraw: self.ax.lines.remove(self.vmaxline[0])
		self.vmaxline=self.ax.plot([self.Vmax.get(),self.Vmax.get()],[ymin,ymax],'--r')
		self.Canvas.draw()

	def OnSetVels(self,event):
		#set UL.Vmin and UL.Vmax based on left/right mouse clicks
		#Also redraw the velocity limit lines on UL.AX
		vel=event.xdata
		ymin,ymax=self.ax.get_ylim()
		if event.button==1:
			self.Vmin.set(vel)
			self.DrawVmin(redraw=True)
		elif event.button==3:
			self.Vmax.set(vel)
			self.DrawVmax(redraw=True)
	def onSaveButton(self):
		#Saves the current values in the GUI to the log file.
		#Determine flag based on true/false values
		flag=1*int(self.OKFlag.get())+2*int(self.BlendFlag.get())+4*int(self.UpperFlag.get())+8*int(self.LowerFlag.get())
		#Save information to the LOG
		self.log[self.zstr,self.ion,self.line,'flag']=str(flag)
		self.log[self.zstr,self.ion,self.line,'vel']=str(self.Vmin.get()),str(self.Vmax.get())
		self.log[self.zstr,self.ion,self.line]=self.Notes.get()
		self.log[self.zstr,self.ion,self.line,'colour']=self.Colour.get()
		#Destroy the widget
		self.Canvas.stop_event_loop()
		self.ULmaster.destroy()
		#Quit UPDATELOG
		return		
def LoadLineList(file):
	""" 
	LOADLINELIST - Loads the supplied linelist into a dictionary for TB_LINEFINDER

	Call - LLIST=LoadLineList(FILE)
	
	INPUT: FILE - The input ASCII file with the linelist for identification.

			Each line in the file must contain information on each spectral line
			of interest in the format (whitespace delimited):
				REST_WAVELENGTH ION SHORT_WL OSCILLATOR_STRENGTH
	
				REST_WAVELENGTH - The rest wavelength of the spectral line (e.g. 1215.6701)
				ION - The species (including ionization state) of the spectral line (e.g. HI)
				SHORT_WL - A short identifier of the wavelength (e.g. 1215)
					Needed to descrimante between absorption lines of a given ION
				OSCILLATOR_STRENGTH - The oscillator strength of the feature
					This is only used for display purposes. (e.g. 0.4164)
	
				EXAMPLE:
				1215.6701 HI 1215 0.4164
	
			If Ly-alpha not present, LOADLINELIST will add it by default
			MUST CONTAIN THE SAME NUMBER OF COLUMNS IN EACH LINE (unless commented)
			Can use '#' at the start of any line to comment

	OUTPUT: LLIST - A dictionary with the linelist information. LLIST is structured as follows:
		LLIST['ions'] - A list of all the ion species in the linelist (ION column)
		LLIST['lines'] - A dictionary (keyed by an ION in LLIST['ions']) of all the
			SHORT_WL identifiers for a given ION
	
		LLIST['lines'][ION] is a list of all the SHORT_WL for a given ION
	
		LLIST[ION,SHORT_WL] - A tuple with (WAVELENGTH,OSCILLATOR_STRENGTH) for 
			the spectral line identified ION,SHORT_WL
	
	NOTES
		For simplicity, SHORT_WL is abbreviated LINE in the code. 
		Will add Ly-alpha (HI 1215) by default if not present in FILE

"""
	#Read in FILE.
	#FILE in format: wvl	ION	lineID	fval	
	data=np.genfromtxt(file,comments='#',dtype=type('str'))
	#Create LLIST, the linelist dictionary
	llist={}
	#Create the list of ions and lines
	llist['ions']=[]
	llist['lines']={}
	#For each line in FILE
	for ii in range(len(data)):
		ion=data[ii][1]
		line=data[ii][2]
		#Convert WAVELENGTH and OSCILLATOR_STRENGTH to a float
		wvl=float(data[ii][0])#WAVELENGTH
		f=float(data[ii][3])#OSCILLATOR_STRENGTH
		#Save tuple to list
		llist[ion,line]=wvl,f
		#Check to see if ION or LINE (aka SHORT_WL) are in LLIST
		#'ions' or 'lines' lists. If not, add them.
		if ion not in llist['ions']:
			llist['ions'].append(ion)
			llist['lines'][ion]=[]
		if line not in llist['lines'][ion]: llist['lines'][ion].append(line)
	#Line list requires Ly-alpha. Add to LLIST if not present.
	if ('HI','1215') not in llist:
		llist['HI','1215']=1215.6701,0.41640
	#REturn LLIST
	return llist	
class LineAdder:
	""" 
	CLASS LINEADDER - The graphical interface for checking all spectral lines in
		linelist (LLIST dictionary) visually, selecting which spectral
		lines to add to the logfile (LOG dictionary). This window contains
		a 3-column list of every available spectral line as a Checkbutton
		for selecting which lines the user wants to add to the LOG dictionary


	Call - LA=LineAdder(RADIOLIST,VELPLOTWIN)

	INPUTS: RADIOLIST - A list of all the spectral lines available to generate
			the Checkbuttons for selection. Each element in RADIOLIST
			is a tuple of (ION,LINE) (the spectrl line key in LLIST dictionary)
		VELPLOTWIN - The TKinter Window that will contain the widget for selecting lines

	ATTRIBUTES:
		LA.ADDLINES - A dictionary with all the spectral lines selected by the Checkbuttons.
			ADDLINES is keyed by every element in RADIOLIST. The value is either TRUE/FALSE
			depending if the corresponding Checkbutton was clicked (TRUE) or not (FALSE)
		LA.MASTER - The master TKinter window (VELPLOTWIN)
		LA.POPUP - The Frame of the TKinter window (i.e. the popup menu)
		LA.VARLIST - A 	dictionary with the variables that contain whether a given
			Checkbutton has been pressed or not. VARLIST is keyed by each element
			in RADIOLIST, and is a TKINTER.BOOLEANVAR.
		LA.SAVEBUTTON - THe widget button for ending the widget to enable saving
			information to LOG dictionary.
		LA.CLOSELAMENU - A function to close VELPLOTWIN and return the class stuff
			back to the place where LINEADDER was called.

	NOTES
		THe only attribute that is really needed is LA.ADDLINES. The rest is just
		for internal use of LINEADDER.
	
		The grid of Checkbuttons will have the same layout as the figure showing
		the velocity plots corresponding to each Checkbutton.
	
		On pressing the SAVEBUTTON, VELPLOTWIN will be destroyed.

	"""
	#Intialize the class when called
	def __init__(self,radiolist,VelPlotWin):
		if debug: print "Initialize LineAdder Function"
		#Define ADDLINES and populate each spectral line key
		#from RADIOLIST with FALSE value (i.e. do not add line)
		self.addlines={}
		for key in radiolist:self.addlines[key]=False
		#Define the parent window, and set up Tkinter frame
		self.master=VelPlotWin
		self.popup = Tkinter.Frame(self.master)
		self.popup.grid()
		self.master.title('LineAdder')
		#Generate list of lines and radio buttons
		self.varlist={}
		#ROW is the row number in the LA.POPUP grid
		row=1
		#NCOL is the number of columns (set to 3) of Checkbuttons
		ncol=3
		#A dummy variable to keep track of which column the next
		#Checkbutton should occupy
		colnum=0
		#FOr each spectral line key, create a Checkbutton
		#and the associated Tkinter variable for that button.
		for key in radiolist:
			self.varlist[key]=Tkinter.BooleanVar()
			self.varlist[key].set(self.addlines[key])
			rb=Tkinter.Checkbutton(self.popup,text='%s (%s)'%key,\
				variable=self.varlist[key])
			rb.grid(column=colnum,row=row,sticky='EW')
			colnum+=1#Increment COLNUM for next Checkbutton
			#If COLNUM is at the maximum number of columns,
			#increment ROW and reset COLNUM
			if colnum==ncol:
				colnum=0
				row+=1
		#If the row does not have enough buttons to fill it, increment
		#ROW
		if colnum>0: row+=1
		#Define the "Save to Log" button for exiting the VELPLOTWIN window
		#This will run the function CLOSELAMENU (below)
		self.savebutton=Tkinter.Button(self.popup,text="Save to Log", command=self.closeLAMenu)#Runs OnFitButton function when clicked
		self.savebutton.grid(column=0,row=row,columnspan=ncol,sticky='EW')
		#Tell VELPLOTWIN to wait until window closes.
		self.popup.wait_window()

	#Function to close the LineAdder menu (i.e. VELPLOTWIN)
	def closeLAMenu(self):
		#Obtain the value of each spectral line Checkbutton
		for key in self.varlist:
			val=self.varlist[key].get()
			#If the value is 1 (i.e. on), update LA.ADDLINES
			#to TRUE for the given spectral line key.
			if val==1:
				self.addlines[key]=True
		if debug: print "LineAdder: Closing Menu"
		#Destroy the TKINTER window VELPLOTWIN and it's frame.
		self.popup.destroy()
		self.master.destroy()
		return
def VelPlots(log,z,llist,fits,VelPlotWin,VPfig,VelFig):
	""" 
	VELPLOTS - The function that plots all lines for a provided redshift 
		within provided linelist LLIST as a velocity profile. It will
		then spawn a Tkinter window to select the spectral lines to
		add to the LOG dictionary using the LINEADDER
		class (defined above). Upon line selection, the user will be queryed
		about the COLOUR,FLAG, VMIN/VMAX and NOTE for each spectral line added
		to the log with the GETUSERINPUT.

	Call - LOG=VelPlots(LOG,Z,LLIST,FITS,VELPLOTWIN,VPFIG, VELFIG)

	INPUTS: LOG - The LOG dictionary defined from READLOG with the
			spectral information
		Z - The redshift of the system of interest (float)
		LLIST - The linelist dictionary LLIST from input linelist file
		FITS - The name of the input spectrum ASCII file to be plotted
		VELPLOTWIN - The TKinter window that contains the matplotlib Canvas
		VELFIG - The matplotlib TKagg canvas for velocity profiles in VELPLOTWIN window
		VPFIG - The matplotlib figure embedded in the VELFIG canvas

	OUTPUT: LOG - The edited LOG dictionary inputted into VELPLOTS.

	NOTES:
		The layout of the VELPLOTWIN canvas should be in the same order
		as the Checkbuttons that appear in the LINEADDER GUI window. It
		will have 3 columns total, and as many rows as needed for each 
		spectral line within the wavelength range of the spectrum.

	"""
	#Load the spectrum from the ascii file with SPECFITS
	#IWVLNGTH and ISPECTRUM are numpy arrays with the
	#wavelength and flux of the spectrum. 
	iwvlngth, ispectrum=specfits(fits)
	if debug: print "VelPlots Loaded spectrum"
	#LKEYS will contain a list of all the spectral
	#line keys within LLIST that are within the 
	#wavelength range of the spectrum.
	lkeys=[]
	#Sort the ions list alphabetically to put
	#IONS in order
	llist['ions'].sort()
	#Loop through each ION/LINE combination in the linelist
	#to identify which lines are covered by the spectrum.
	#If that is true, add the LLIST key to LKEYS.
	for ion in llist['ions']:
		#Sort the lines for a given ION so they are in order
		llist['lines'][ion].sort()
		for line in llist['lines'][ion]:
			#DOuble check ION/LINE  is in LLIST
			#If this fails, something is fundamentally flawed...
			if (ion,line) in llist:
				#Get the wavelength of the spectral line
				wvl,f=llist[ion,line]
				#Redshift the line....
				wvl=wvl*(1.0+z)
				if debug: print "Wavelength check: is %.3f within (%.3f,%.3f)"%(wvl,min(iwvlngth),max(iwvlngth))
				#If redshifted line covered by spectrum, add to LKEYS
				if wvl>min(iwvlngth) and wvl<max(iwvlngth):
					if debug: print "Passed wavelngth check",ion,line
			 		lkeys.append((ion,line))
	#NCOL is the number of columns in the VELPLOTWIN velocity profile grid
	ncol=3
	#NROW is the number of rows in the VELPLOTWIN grid.
	nrow=len(lkeys)/ncol+1
	if debug: print "VelPlots panel Nrow,ncol,nkeys:", nrow,ncol,len(lkeys)
	#Set the limits of the plotted velocity profile to be +/- 1000 km/s
	vmin=-1000
	vmax=1000
	if debug: print "VelPlots Plotting lines:",lkeys
	#RADIOLIST is the list of LLIST spectral lines that require
	#Checkbuttons in the LINEADDER GUI window
	radiolist=[]
	#The axes instance variable to share with all other panels
	shareax=None
	#Loop through all spectral lines for plotting and plot them!
	for ii in range(len(lkeys)):
		#Get the spectral line key (ION,LINE) tuple
		lkey=lkeys[ii]
		#Grabd the wavelenght and oscilaltor strength of the line
		wvl,f=llist[lkey]
		#WVL should be Redshifted 
		wvl=wvl*(1.0+z)
		#IND refers to the matplotlib subplot index for the
		#velocity profile
		ind=ii+1
		#Get the spectrum slice within vmin/vmax of the spectrum
		#using VELSPEC (centred at the redshifted wavelength WVL)
		#TMPVEL is the velocity array, TMPSPEC is the corresponding flux
		tmpvel,tmpspec=velspec(wvl, iwvlngth, vmin,vmax, ispectrum)
		#Variable to contain axes instance for the subplot
		ax=None
		#To share the same zoom on all plotting windows, the first window
		#must be identified (ii==0), while subsequent widnows must share
		#the x/y axes (SHAREX/SHAREY)
		if ii==0:
			#Add a subplot to the VELPLOTWIN figure, call it AX
			ax=VPfig.add_subplot(nrow,ncol,ind)
			shareax=ax
		else:
			ax=VPfig.add_subplot(nrow,ncol,ind,sharex=shareax)
		
		#Plot the velocity profile to AX
		ax.plot(tmpvel,tmpspec,'k',drawstyle='steps')
		#Add key to RADIOLIST
		radiolist.append(lkey)
		#GEt y-axis limits of the velocity profile, and modify
		#YMIN to include and extra 10% of the original heaigh
		ymin,ymax=ax.get_ylim()
		ymin=ymin-0.1*(ymax-ymin)
		#get the ION,LINE identifier from the spectral line key
		ion,line=lkey
		#Set the title of the subplot to the idetifier, and 
		#provide the oscillator strength for reference
		title='%s %s\nf=%.5f'%(ion,line,f)
		ax.set_title(title)
		#Loop through the current LOG dictionary, and plot
		#a vertical dashed line to inform the user if there
		#is potentially another line from a different system.
		#LZ is the redshift of a system in the LOG allready
		for lz in log['zs']:
			#LION is the ion in the log
			for lion in log['ions']:
				#LLINE is the line identifier for LION
				for lline in log['lines'][lion]:
					#If the log has the key (LZ,LION,LLINE) and (LION,LLINE) are in the linelist LLIST
					if (lion,lline) in llist and (lz,lion,lline) in log:
						#Get the wavelength of the LION,LLINE
						#And redshift it by LZ
						wl,f=llist[lion,lline]
						#Get the velocity of this line
						#with respect to the the spectral
						#Line being plotted.
						vel=(wl*(1.0+float(lz))-wvl)/wvl*2.998E5
						#If it's ALSO within the velocity
						#profile of the spectral line plotted
						if vel>vmin and vel<vmax:
							#Plot a vertical dashed line, and add
							#a label to inform the user of the contaminating
							#line.
							if debug: 'VelPlot Adding Log redshift',lz
							col=log[lz,lion,lline,'colour']
							ax.plot([vel,vel],[ymin,ymax],'--'+col, linewidth=3)
							fval='%.5f'%f
							label='z=%s\n%s %s'%(lz,lion,lline)
							ax.text(vel,ymax,label,color=col,va='top',ha='left',rotation='vertical')
		#Set the subplot velocity limits to the velocity range specified by VMIN/VMAX
		ax.set_xlim(vmin,vmax)
	#Update the VELFIG canvas
	VelFig.draw()
	if debug: print "VelPlots RadioList:", radiolist
	#Use the LINEADDER class to generate the GUI window to select lines
	#for adding to LOG
	lineadder=LineAdder(radiolist,VelPlotWin)
	#ADDLINES is a dicitionary keyed by the elements of RADIOLIST.
	#Each key corresponds to whether or not the line should be added
	#(true if add, false if not)
	addlines=lineadder.addlines
	if debug: print "VelPlots addlines:", addlines
	#COnvert the redshift to a string for the log. Should be 5 decimcal places.
	zstr='%.5f'%z
	#Loop through each spectral line  within ADDLINES
	for key in addlines:
		#If the user selected it (i.e. value==True)
		if addlines[key]==True:
			#Check if that system allready exists. If not
			#add ZSTR to the list of redshifts in LOG['zs']
			if zstr not in log['zs']: log['zs'].append(zstr)
			#Obtain the ION,LINE identifier of the spectral line
			ion,line=key
			#Check to see if ION allready in LOG['ions'] list.
			#If not, add it, and make an entry into LOG['lines'] dictionary
			#for that ion.
			if ion not in log['ions']:
				log['ions'].append(ion)
				log['lines'][ion]=[]
			#IF LINE identifier not in the LOG['lines'][ION] list, add it
			if line not in log['lines'][ion]: log['lines'][ion].append(line)
			#Pass the system/spectral line information to GETUSERINPUT
			#to query user for log information about the line, and update
			#the LOG dictionary appropriately
			wvl,f=llist[ion,line]
			UL=UpdateLog(log, z,ion,line,wvl, iwvlngth, ispectrum)
			log=UL.log
	#Now you have an updated LOG dictionary for a given system redshift, return it.
	return log

class linefinder_tk(Tkinter.Tk):# Base for standard window
	""" 
	CLASS LINEFINDER_TK is the main TKINTER interface for TB_LINEFINDER. IT grabs the input/output file names
		and initializes the program by plotting the spectrum. Once loaded, the user can use keyboard
		shortcuts to identify spectral lines (from the inputted linelist) to a logfile
		which is sotred in the LOG dictionary

	CALL  LF=linefinder_tk(None)

	ATTRIBUTES
		LF.root - The Tkinter window identifier drop menu, and buttons
		LF.entryIn - TKINTER string vairable that contains the input spectrum file name
		LF.fits - The spectrum file name (string) that is obtained from LF.entryIn
		LF.entryOut - TKINTER string variable that contains the output log file name	
		LF.logfile - The logfile filename string obtained from LF.entryOut
		LF.entryLlist - TKINTER string variable with input linelist filename
		LF.log - the LLOG dictionary that contains the data to be written to LF.entryLlist
		LF.SProot - The TKInter parent for the full spectrum plot window
		LF.SpecPlot - The matplotlib canvas embedded in LF.SPROOT
		LF.ax - the Matplotlib Axes instance for the full spectrum plot window (LF.SpecPlot)
		LF.fig - Matplotlib figure instance for the full spectrum plot window (LF.SpecPlot)

		LF.onExit() - function to close the LF window and exit TB_LINEFINDER for good
		LF.initialize() - The initilization functuon of the LF window (set up the user input fields,
		LF.PlotSpec() - Function to plot the spectrum from LF.fits
		LF.PlotFits() - Function to generate the the TKinter LF.SProot window
		LF.UpdatePlot() - Refresh the LF.SpecPlot canvas with current LOG dictionary content
		LF.LogMenu(LOGMENUWIN) - Defines the TKinter window for editing the current LOG dictionary content
			in the LOGMENUWIN child window.
		LF.OnChangeLog() - Update LOG dictionary content with user input (done via command line)
		LF.OnRemoveLine() - Remove a spectral line from the LOG dictrionary for a certain system
		LF.OnRemoveSystem() - Remove an entire absorption system from the LOG dictionary
		LF.OnExitEdit() - Closes the LOG Editing Menu GUI started in LF.LogMenu()
		LF.UpdateIonMenu() - Update the drop-down Ion menu in the LF.LogMenu to contain only ions in LOG 
		LF.UpdateLineMenu() - Update the drop-down line menu to include lines for a given ion in the LOG
		LF.onQ() - Function when 'Q' is pressed in LF.SpecPlot window (Quits LineFinder)
		LF.onVelPlots(Z) - Sets up the appropriate TKinter windows and runs the VELPLOTS function for a given redshift
		LF.onL(EVENT) - Function when 'L' is pressed in LF.SpecPlot window (Runs LF.onVelPlots for cursour position on Ly-a line)
			LF.onZ() - Function when 'Z' is pressed in LF.SpecPlot window (gets redshift from user in command line and runs LF.onVelPlots)
		LF.onA(EVENT) - Function when 'A' is pressed in LF.SpecPlot window (Runs LF.onVelPlots for cursour position on a given line from command line)
		LF.onE() - Function when 'E' is pressed in LF.SpecPlot window (Runs LF.OnChangeLog)
		LF.onH() - Function when 'H' is pressed in LF.SpecPlot window (displays key options)
		LF.onKey(EVENT) - Obtain LF.SpecPlot event and decide which function to run.
		LF.OnFindLines() - The initialization of the Line Finder application. Sets up LF.SpecPlot and starts listening for interactive events
		LF.getSpecFile() - File browser window to select the input spectrum ascii file
		LF.getLlistFile() - File browser window to select the input linelist file
		LF.getOutFile() - File browser window to select the name of the output log file

	NOTES
		This is the main driver of the GUI interface, and is where new features should be added.
	
		KNOWN BUG - LF.LineMenu() will always complain about missing keys in dictionary on first use of it. Ignore unless repeated.
			The message reads:

			Traceback (most recent call last):
			  File "/usr/lib64/python2.7/lib-tk/Tkinter.py", line 1470, in __call__
			    return self.func(*args)
			  File "tb_linefinder.py", line 1103, in UpdateIonMenu
			    self.ion.set(zions[0])
			  File "/usr/lib64/python2.7/lib-tk/Tkinter.py", line 1826, in __getattr__
			    return getattr(self.tk, attr)
			AttributeError: ion
			Exception in Tkinter callback
			Traceback (most recent call last):
			  File "/usr/lib64/python2.7/lib-tk/Tkinter.py", line 1470, in __call__
			    return self.func(*args)
			  File "tb_linefinder.py", line 1114, in UpdateLineMenu
			    for line in self.log['lines'][self.ion.get()]:
	
	
	"""
        #This defines the main window
        def __init__(self,parent):
		#Root identifier for Tkinter window
                self.root=Tkinter.Tk.__init__(self,parent)
                #self.parent =parent
                self.initialize()
	#Function to quit the entire GUI
        def onExit(self):
                print "Quitting..."
                self.quit()
	def onTutorial(self):
		global usetutorial
		if usetutorial:
			usetutorial=False
			tkMessageBox.showinfo("Help Message", "Tutorial mode is now off.")
		else:
			usetutorial=True
			tkMessageBox.showinfo("Help Message", "Tutorial mode is now on.")
		return
	#Browser dialog buttons for input/output files...
	#... To get the Input spectrum file
	def getSpecFile(self):
		#Options for the browser dialog window
		browser_opt={}
		#Have access to all files
		browser_opt['filetypes']=[('all files','.*')]
		#Set initial file to current spectrum file
		browser_opt['initialfile']=self.specIn.get()
		#Label the window
		browser_opt['title']='Select Input spectrum (ascii)'
		#Run the dialog, and save the output to the spectrum file name
		self.specIn.set(tkFileDialog.askopenfilename(**browser_opt))
		return
	#... To get the Input linelist file
	def getLlistFile(self):
		#Options for the browser dialog window
		browser_opt={}
		#Have access to all files
		browser_opt['filetypes']=[('all files','.*')]
		#Set initial file to current linelist file
		browser_opt['initialfile']=self.entryLlist.get()
		#Label the window
		browser_opt['title']='Select Line List'
		#Run the dialog, and save the output to the llist file name
		self.entryLlist.set(tkFileDialog.askopenfilename(**browser_opt))
		return
	#... To get the output logfile
	def getOutFile(self):
		#Options for browser dialog window
		browser_opt={}
		#Have access to all files
		browser_opt['filetypes']=[('all files','.*')]
		#Label the window
		browser_opt['title']='Name of output file'
		#Set initial file to current ouput file
		browser_opt['initialfile']=self.entryOut.get()
		#Run the dialog, and save the output to the logfilename
		self.entryOut.set(tkFileDialog.askopenfilename(**browser_opt))
		return


	#Function for initilizaing the main GUI window
        def initialize(self):
		#Open tutorial message box saying what to do
		if usetutorial: tkMessageBox.showinfo("Help Message", "Select input spectrum and linelist file, and name output logifle.")
                self.grid()
                #Drop down menu for main GUI (e.g. exit button)
                mb=Tkinter.Menubutton(self,text='Menu',relief="raised")
                mb.grid(column=0,row=0, sticky='EW')
                picks=Tkinter.Menu(mb,tearoff=0)
                picks.add_command(label="Tutorial mode on/off",command=self.onTutorial)
                picks.add_command(label="Exit",command=self.onExit)
                mb.config(menu=picks)

                #########################
                #Define the Input fields#
                #########################

                #The label for the spectrum file
                labelIn=Tkinter.StringVar()#Name of string in label
                labelin=Tkinter.Label(self,textvariable=labelIn,\
                        anchor="w",fg="black")
                labelin.grid(column=0, row=1, sticky='EW')
                labelIn.set(u"Input spectrum:")#Initial value of Variable

                #The field for the spectrum file
                self.specIn=Tkinter.StringVar()
                entryin=Tkinter.Entry(self,textvariable=self.specIn)#,anchor="w", fg="black")
                entryin.grid(column=1,row=1,stick='EW')
                self.specIn.set(initinspec)
		
		#Button for input spectrum file browser
		specbutton=Tkinter.Button(self,text="Search",command=self.getSpecFile)
		specbutton.grid(column=2,row=1, sticky='EW')

                #The label for the linelist
                labelLlist=Tkinter.StringVar()#Name of string in label
                labelllist=Tkinter.Label(self,textvariable=labelLlist,\
                        anchor="w",fg="black")
                labelllist.grid(column=0, row=2, sticky='EW')
                labelLlist.set(u"Linelist:")#Initial value of Variable

                #The field for the linelist
                self.entryLlist=Tkinter.StringVar()
                entryllist=Tkinter.Entry(self,textvariable=self.entryLlist)#,anchor="w", fg="black")
                entryllist.grid(column=1,row=2,stick='EW')
                self.entryLlist.set(initllist)

		#Button for input linelist browser
		llistbutton=Tkinter.Button(self,text="Search",command=self.getLlistFile)
		llistbutton.grid(column=2,row=2, sticky='EW')
		

                #The label for the output logfile
                labelOut=Tkinter.StringVar()#Name of string in label
                labelout=Tkinter.Label(self,textvariable=labelOut,\
                        anchor="w",fg="black")
                labelout.grid(column=0, row=3, sticky='EW')
                labelOut.set(u"Output logfile:")#Initial value of Variable

                #The field for the output logfile
                self.entryOut=Tkinter.StringVar()
                entryout=Tkinter.Entry(self,textvariable=self.entryOut)#,anchor="w", fg="black")
                entryout.grid(column=1,row=3,stick='EW')
                self.entryOut.set(initlog)

		#Button for output logfile browser
		outlogbutton=Tkinter.Button(self,text="Search",command=self.getOutFile)
		outlogbutton.grid(column=2,row=3, sticky='EW')


                #Define button for Finding lines
                plotbutton=Tkinter.Button(self,text=u'Find Lines!',\
                        command=self.OnFindLines)#Runs OnFitButton function when clicked
                plotbutton.grid(column=0,row=4,columnspan=3,sticky='EW')

		#Initialize the LOG dictionary
		self.log=None
	#Function to Plot the input spectrum
	def PlotSpec(self):
		#Load spectrum and plot it to SpecPlot window
		iwvlngth, ispectrum=specfits(self.fits)
		self.ax.plot(iwvlngth,ispectrum,'k',drawstyle='steps')
		self.SpecPlot.draw()
		#The default scaling removes the +/- 2.5 percentile outliers.
		#This number was set arbitrarily for practice spectrum, and might not be appropriate
		ymin=np.percentile(ispectrum,2.5)
		ymax=np.percentile(ispectrum,97.5)
		#Set the y axis limits based on removing outliers
		self.ax.set_ylim(ymin,ymax)
		if debug: print "PlotSpec flux limits: ", ymin, ymax
		return
	#Function to load the Tkinter window for plotting the full spectrum
	def PlotFits(self):
		#Generate the figure for plotting
		self.fig=plt.figure()
		#Generate the Tkinter window
		self.SProot=Tkinter.Tk()
		#Name of TKInter window
		self.SProot.wm_title("Full Spectrum")
		#Create a TkAgg canvas for plotting using self.fig
		self.SpecPlot=FigureCanvasTkAgg(self.fig,master=self.SProot)
		self.SpecPlot.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
		#Add a toolbar for zooming, scrolling, resetting, etc of the plot window
		NavSpecPlot=NavigationToolbar2TkAgg(self.SpecPlot, self.SProot)
		#Generate an matplotlib axis for plotting the spectrum
		self.ax=plt.subplot(1,1,1)
		#Plot the spectrum
		self.PlotSpec()
		return


	#The function related to updating the full specturm plot with the log information
	#This includes coloured vertical lines for each identified feature and the label
	#the redshift and absorption feature at that line
	def UpdatePlot(self):
		#Clear the MPL axis and reset the canvas
		self.ax.clear()
		self.SpecPlot.draw()
		self.fig.set_canvas(self.SpecPlot)
		#Replot the spectrum from input file
		self.PlotSpec()
		if debug: print "Updating Plot", self.log
		#Get the y-limits of the plot for making veritcal lines
		ymin,ymax=self.ax.get_ylim()
		#Go through each system and plot each absorption feature in the log.
		#Use the colour identified in LOG[Z,ION,LINE,'colour']
		#Label with Z, ION, LINE
		for z in self.log['zs']:
			if debug: print "UpdatePlot redshift", z
			for ion in self.log['ions']:
				if debug: print "\tUpdatePlot ion", ion
				for line in self.log['lines'][ion]:
					if debug: print "\t\tUpdatePlot line", line
					if debug: print "\t\t In linelist:", (ion,line) in self.llist
					if debug: print "\t\t In log:", (z,ion,line) in self.log
					if (ion,line) in self.llist and (z,ion,line) in self.log:
						if debug: print "\t\t\tUpdatePlot", z, ion, line
						wl,f=self.llist[ion,line]
						wl=wl*(1.0+float(z))
						#By default use black if no colour in LOG
						col='k'
						if (z,ion,line,'colour') in self.log:
							col=self.log[z,ion,line,'colour']
						self.ax.plot([wl,wl],[ymin,ymax],':'+col)
						fval='%.5f'%f#Oscillator strength string if necessary
						label='z=%s\n%s %s'%(z,ion,line)
						self.ax.text(wl,ymax,label,color=col,va='top',ha='left',rotation='vertical')
					#else: print "Warning: Line %s %s not in linelist."%(ion,line)
		#Redraw the canvas to update the plotting window
		self.SpecPlot.draw()
		if debug: print "UpdatePlot Draw Figure" 
		return
	# Get user input about updating the LOG file in LOGMENU using the command line
	def OnChangeLog(self):
		if debug: print "Running OnChangeLog"
		#Grab the redshift (ZZ), ion (IION), and line identifier (LLINE) in the LOG keys for updating
		zz=self.z.get()
		iion=self.ion.get()
		lline=self.line.get()
		wvl,f=self.llist[iion,lline]
		#Grab spectrum, and send to UpdateLog window for log editing
		iwvlngth, ispectrum=specfits(self.fits)
		UL=UpdateLog(self.log, float(zz),iion,lline,wvl, iwvlngth, ispectrum)
		#Save the updated log file
		self.log=UL.log
		return
	#What to do when the "Remove system" button is pressed in the LogMenu widget
	def OnRemoveSystem(self):
		if debug: print "Running OnRemoveSystem"
		#Obtain the redshift of the desired removed system selected via the drop menu 
		rmz=self.z.get()
		#go through each ION/LINE combination, and see if it has an entry for the system at redshfit RMZ 
		for rmion in self.log['ions']:
			for rmline in self.log['lines'][rmion]:
				#If combination is in the LOG, remove it
				if (rmz,rmion,rmline) in self.log:
					self.log.pop((rmz,rmion,rmline),0)
					self.log.pop((rmz,rmion,rmline,'flag'),0)
					self.log.pop((rmz,rmion,rmline,'colour'),0)
					self.log.pop((rmz,rmion,rmline,'vel'),0)
		#Figure out which index in the list of redshifts the RMZ absorber corresponds to,
		#and remove from the redshift list
		inds=np.where(rmz==self.log['zs'])[0]
		for ind in inds:
			self.log['zs'].pop(ind)
		#Check to see if there are any ions/lines that are no longer present in the LOG.
		#if so, they need to be removed from the list of possible lines and ions.
		for rmion in self.log['ions']:
			removeion=True#Boolean to see if whether or not the ion needs to be removed
			for rmline in self.log['lines'][rmion]:
				removeline=True#boolean to see if the line of the ion needs to be removed
				for rmz in self.log['zs']:
					#if the ion/line combination is in LOG all ready after removing redshift, we need to keep it.
					if (rmz,rmion,rmline) in self.log:
						#Set the boolean values to false so ION/line combo is not removed
						removeline=False
						removeion=False
				#Get rid of line from list if deemed apprpriate to remove
				if removeline:
					inds=np.where(rmline==self.log['lines'][rmion])[0]
					for ind in inds:
						self.log['lines'][rmion].pop(ind)
			#If the list of possible lines is empty and we have said we should remove the ion,
			#remove the ion appropriately
			if len(self.log['lines'][rmion])<1 and removeion:
				inds=np.where(rmion==self.log['ions'])[0]
				for ind in inds:
					self.log['ions'].pop(ind)
		#Update the drop menus appropriately to encompass any changes
		self.UpdateLineMenu()
		self.UpdateIonMenu()
		self.UpdateRedshiftMenu()
		return
	#if a specific line for a given absorption system is desire to be removed, get rid of the
	#specific entry in the LOG (including flags, colours, and velocity information)
	def OnRemoveLine(self):
		if debug: print "Running OnRemoveLine"
		#Grab the current redshift, ion, and line ID for removal
		rmz=self.z.get()
		rmion=self.ion.get()
		rmline=self.line.get()
		#Double check the corresponding key is in the dictionary, if so remove it
		if (rmz,rmion, rmline) in self.log:
			self.log.pop((rmz,rmion,rmline),0)
			self.log.pop((rmz,rmion,rmline,'flag'),0)
			self.log.pop((rmz,rmion,rmline,'colour'),0)
			self.log.pop((rmz,rmion,rmline,'vel'),0)
		#Update the drop-down menus
		self.UpdateIonMenu()
		self.UpdateLineMenu()
		return
	#Save the log file and exit the LogMenu editing widget
	def OnExitEdit(self):
		WriteLog(self.log,self.logfile,self.llistfile,self.fits)
		#Close the widgets
		self.EditPopup.destroy()
		self.EditMaster.destroy()
		return
	#Update the drop-down menu displayin the ions in the LogMenu list
	def UpdateIonMenu(self, *args):
		zions=[]
		#Comb LOG dictionary for every ion in LOG, and add it to a list of possible ions
		for ion in self.log['ions']:
			for line in self.log['lines'][ion]:
				if ion not in zions:
					if (self.z.get(),ion,line) in self.log:
						zions.append(ion)
		if debug: print "UpdateIonMenu set self.ion",zions[0]
		#Set the default value to <ion>
		self.ion.set("<ion>")
		#update the IonMenu by removing all previous entries and adding the new
		menu=self.IonMenu['menu']
		menu.delete(0,'end')
		for ion in zions:
			menu.add_command(label=ion, command=lambda ion=ion: self.ion.set(ion))
		return
	#Update drop-down menu displaying the lines for the ION selected in the Ion drop-down menu in LogMenu
	def UpdateLineMenu(self, *args):
		zlines=[]
		#Comb LOG dictionary for every line for the given ion in the LOG and add to list of possible lines
		if self.ion.get() in self.log['lines']:
			for line in self.log['lines'][self.ion.get()]:
				if (self.z.get(),self.ion.get(),line) in self.log:
					zlines.append(line)
		if debug: print "UpdateLineMenu linelist",zlines
		#Set the default value to <line>
		self.line.set('<line>')
		#Clear drop down meny and add new line list to it.
		menu=self.LineMenu['menu']
		menu.delete(0,'end')
		for line in zlines:
			menu.add_command(label=line, command=lambda line=line: self.line.set(line))
		return


	#Update drop-down menu displaying the redshifts in the drop-down menu in LogMenu
	def UpdateRedshiftMenu(self, *args):
		zs=[]
		#Comb LOG dictionary for all redshifts
		for z in self.log['zs']:
			zs.append(z)
		if debug: print "UpdateRedshiftMenu list"
		#Set the default value to first redshift
		self.z.set('<z>')
		#Clear drop down meny and add new redshifts to list
		menu=self.ZMenu['menu']
		menu.delete(0,'end')
		for z in zs:
			menu.add_command(label=z, command=lambda z=z: self.z.set(z))

		return


	#Function setting up the log-editing menu, with a list of absorbtion redshifts, ions, and lines.
	#Options include editing the log of a given entry, remove a system redshift, or a line at a given redshift.
	def LogMenu(self,LogMenuWin):# Base for Log Menu
	        #This defines the Log Menu window
		self.EditMaster=LogMenuWin
		self.EditMaster.title('Log Editing Menu')
		self.EditPopup = Tkinter.Frame(self.EditMaster)
		self.EditPopup.grid()
		if debug: print "Loading Log editing Menu"
		#System redshift dropdown menu
		labelzvar=Tkinter.StringVar()
		labelzvar.set('z')
		labelz=Tkinter.Label(self.EditPopup,textvariable=labelzvar,anchor="w",fg="black")
		labelz.grid(column=0,row=1,sticky='EW')
		self.z= Tkinter.StringVar(self.EditPopup)
		#For the selected redshift, update the dropdown ion meny accordingly
		self.z.trace('w',self.UpdateIonMenu)
		#Sort redshifts numerically (starting with the lowest)
		self.log['zs'].sort()
		self.z.set('<z>') # default value
		self.ZMenu = apply(Tkinter.OptionMenu, (self.EditPopup, self.z) +tuple(self.log['zs']))
		self.ZMenu.grid(column=1,row=1,stick='EW')
		#Example code for changing menus (see sample in UpdateLineMenu for other half); don't use.
		""" 
		    def __init__(self, master):
		        tk.Frame.__init__(self, master)
		        self.dict = {'Asia': ['Japan', 'China', 'Malasia'],
	                     'Europe': ['Germany', 'France', 'Switzerland']}
		        self.variable_a = tk.StringVar(self)
		        self.variable_b = tk.StringVar(self)
		        self.variable_a.trace('w', self.updateoptions)
		        self.optionmenu_a = tk.OptionMenu(self, self.variable_a, *self.dict.keys())
		        self.optionmenu_b = tk.OptionMenu(self, self.variable_b, '')
		        self.variable_a.set('Asia')

		"""


		#Ion dropdown menu
		labelionvar=Tkinter.StringVar()
		labelionvar.set('Ion')
		labelion=Tkinter.Label(self.EditPopup,textvariable=labelionvar,anchor="w",fg="black")
		labelion.grid(column=0,row=2,sticky='EW')
		self.ion= Tkinter.StringVar(self.EditPopup)
		#Need to update the Line Menu
		self.ion.trace('w',self.UpdateLineMenu)
		self.ion.set('') # default value
		self.IonMenu = Tkinter.OptionMenu(self.EditPopup,self.ion,'')
		self.IonMenu.grid(column=1,row=2,stick='EW')

		#Line dropdown menu
		labellvar=Tkinter.StringVar()
		labellvar.set('Line')
		labelline=Tkinter.Label(self.EditPopup,textvariable=labellvar,anchor="w",fg="black")
		labelline.grid(column=0,row=3,sticky='EW')
		self.line= Tkinter.StringVar(self.EditPopup)
		self.line.set('') # default value
		self.LineMenu = Tkinter.OptionMenu(self.EditPopup,self.line,'')
		self.LineMenu.grid(column=1,row=3,stick='EW')

		#Buttons for editing LOG, removing line, removing system, and quitting edit log menu
		#Each should run the appropriate function:
		#	Edit logfile - OnChangeLog
		Editbutton=Tkinter.Button(self.EditPopup,text="Edit Logfile", command=self.OnChangeLog)
		Editbutton.grid(column=0,row=4,columnspan=2,sticky='EW')
		#	Remove Line - OnRemoveLine
		Removebutton=Tkinter.Button(self.EditPopup,text="Remove Line", command=self.OnRemoveLine)
		Removebutton.grid(column=0,row=5,columnspan=2,sticky='EW')
		#	Remove System - OnRemoveSystem
		Systembutton=Tkinter.Button(self.EditPopup,text="Remove System", command=self.OnRemoveSystem)
		Systembutton.grid(column=0,row=6,columnspan=2,sticky='EW')
		#	Quit - OnExitEdit
		Quitbutton=Tkinter.Button(self.EditPopup,text="Quit", command=self.OnExitEdit)
		Quitbutton.grid(column=0,row=7,columnspan=2,sticky='EW')
		#Wait for closing the window before proceeding.
		self.EditMaster.wait_window()

				
	#Function to quit main LineFinder routine.
	def onQ(self):
		#Save log file
		WriteLog(self.log,self.logfile,self.llistfile,self.fits)
		#Close event loop
		self.SpecPlot.stop_event_loop()
		#Destroy widget
		self.SProot.destroy()
		print "Quitting Line Finder"
		return
	#Function to run when you need velocity profiles for a given redshift.
	def onVelPlots(self,z):
		#Open tutorial message box saying what to do
		if usetutorial:tkMessageBox.showinfo("Help Window", "Using the velocity profiles, select which lines are associated with the system.")
		#display velocity plots for given line
		if debug: print "Running VELPLOTS", self.log, z, self.llist, self.fits
		#Create a matplotlib figure
		VPfig=plt.figure(figsize=(8,8))
		#BIGAX is the full axes instance for the entire figure.
		bigax = VPfig.add_subplot(111)    # The big subplot
		#BIGAX is used to set up one label for the x and y axis....
		#All the tick marks, etc, need to be turned off
		bigax.spines['top'].set_color('none')
		bigax.spines['bottom'].set_color('none')
		bigax.spines['left'].set_color('none')
		bigax.spines['right'].set_color('none')
		#bigax.xaxis.tick_labels([])
		#bigax.yaxis.tick_labels([])
		bigax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
		#Label axes in BIGAX
		bigax.set_xlabel(r'Relative Velocity (km s$^{-1}$)')
		bigax.set_ylabel(r'Flux')
		#Set up the spacing in between subplots
		VPfig.subplots_adjust(wspace=0.3, hspace=0.5)
		#Make a widget for selecting velocity profiles to include in LOG
		VelPlotWin=Tkinter.Toplevel(self.root)
		#Create a new widget for displaying velocity profiles with matplotlib canvas
		VProot=Tkinter.Tk()
		VProot.wm_title("Velocity Profiles")
		VelFig=FigureCanvasTkAgg(VPfig,master=VProot)
		VelFig.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
		#Add a matplotlib toolbar for zooming, panning, etc
		NavVelFig=NavigationToolbar2TkAgg(VelFig, VProot)
		VPfig.set_canvas(VelFig)
		if debug: "onVelPlots VPfig figure:", VPfig
		#Run the VELPLOTS function
		self.log=VelPlots(self.log,z,self.llist,self.fits,VelPlotWin,VPfig,VelFig)
		#Save notes to log
		WriteLog(self.log,self.logfile,self.llistfile,self.fits)
		VProot.destroy()
		return
	#Function for what to do when 'l' is pressed on keyboard
	#Set up to assume that location of keystroke is a Ly-alpha line
	#Will find redshift and run OnVelPlots.
	def onL(self,event):
		#Get wavelength from where L was pressed on Full Spectrum window
		wvl=event.xdata
		#look up the rest wavelength of Ly-a
		wl,f=self.llist['HI','1215']
		#Caclulate the redshift
		z=(wvl/wl)-1
		if debug: print "Pressed L, running onVelPlots"
		#Run onVelPlots and save to log
		self.onVelPlots(z)
		#WriteLog(self.log,self.logfile,self.llistfile,self.fits)
		return
	#Function to run when 'z' ia pressed.
	#This will ask user for a redshift, see if it is a float
	#if so, pass the redshift to OnVelPlots.
	def onZ(self):
		z=float(raw_input("What redshift do you need: "))
		if IsFloat(z):	self.onVelPlots(z)
		else: print "Invalid Redshift, try again"
		#WriteLog(self.log,self.logfile,self.llistfile,self.fits)
		return
	#Function when 'a' is pressed.
	#This will ask the user for a given ion/line ID to calculate a redshift
	#From cursor position.
	def onA(self,event):
		#Open tutorial message box saying what to do
		if usetutorial:tkMessageBox.showinfo("Help Window", "Use the terminal window to select which line to use.")
		#Get wavelength position of cursor on Full Spectrum window
		wvl=float(event.xdata)
		ion=None#ION name placeholder
		line=None#LINE ID name placeholder
		badion=True#Boolean to see if the ion selected is in linelist(TRUE)
		badline=True#boolean to see if the line ID selecte is in linelist(TRUE)

		#Ask user for a valid ION from the linelist (complete list is displayed)
		#Will exit once a proper ION is submitted
		while badion:
			print "IONS: ",self.llist['ions']
			ion=raw_input("What ION do you need: ")
			if ion in self.llist['ions']:
				badion=False
			else: print "Bad ION, try again"
		#Do the same, but look for a valid line ID based on ion chose.
		while badline:
			print "LINES: ",self.llist['lines'][ion]
			line=raw_input("What LINE do you need: ")
			if line in self.llist['lines'][ion]:
				badline=False
			else: print "Bad LINE, try again"
		#Check to see if in linelist (SOMETHING IS WRONG IF NOT)
		if (ion,line) in self.llist:
			#Grab the restwavelength, and get redshift
			wl,f=self.llist[ion,line]
			z=(wvl/wl)-1
			if debug: print "OnA redshift",z
			#Pass redshift onto onVelPlots
			self.onVelPlots(z)
		else: print "ION,LINE not in Line List. Try again"
		#WriteLog(self.log,self.logfile,self.llistfile,self.fits)
		return
	#What to do if 'e' is pressed
	#This will open the Log editing menu
	def onE(self):
		#Open tutorial message box saying what to do
		if usetutorial:tkMessageBox.showinfo("Help Window", "Use the dropdown menu to edit logfile.")
		#Create log editing menu child widget
		LogMenuWin=Tkinter.Toplevel(self.root)
		#Send this to the Logmenu routine
		self.LogMenu(LogMenuWin)
		#Close widget when done
		LogMenuWin.destroy()
		#Save information to log
		WriteLog(self.log,self.logfile,self.llistfile,self.fits)
		return
	#What to do if 'h' is pressed'
	#This will display a list of keys to press in terminal
	def onH(self):
		print 'q - quit loop and save logfile'
		print 'a - add system line (on cursor position and command line entry)'
		print 'l - select lya line (on cursor position & show potential absorption)'
		print 'z - display potential absorption for input redshift (command line entry)'
		print 'e - edit log, colours, notes (interface)'
		print 'h - help'
		return
	#Function to figure out, based on a keyboard event, what function to run
	#If key the key will change the log file, it needs to run the Full Spectrum plot update
	def onKey(self,event):
		if debug: print "Event: ", event.key, event.xdata,event.ydata
		if event.key=='q': self.onQ()
		elif event.key=='a': self.onA(event)
		elif event.key=='l': self.onL(event)
		elif event.key=='z': self.onZ()
		elif event.key=='e': self.onE()
		elif event.key=='h': self.onH()
		#Update spectrum window if something has changed
		if event.key not in ['q','h']: self.UpdatePlot()
		return
	#What to do when the Find Lines button is pressed on the main widget
	def OnFindLines(self):
		#Load log file (if it exists)
		self.logfile=self.entryOut.get()
		self.log=LoadLog(self.logfile)
		#Load linelist
		self.llistfile=self.entryLlist.get()
		self.llist=LoadLineList(self.llistfile)
		if debug: print "LLIST:",self.llist
		#Open Figure
		self.fits=self.specIn.get()
		self.PlotFits()
		self.UpdatePlot()
		if debug: print "FIGURE:",self.SpecPlot
		if debug: print "Displaying Figure"
		#For help, display potential keys to start
		self.onH()
		#Show plotting window and wait for keys to be pressed indefinetly
		self.SpecPlot.show()
		self.SpecPlot.mpl_connect('key_press_event',self.onKey)
		#Open tutorial message box saying what to do
		if usetutorial:tkMessageBox.showinfo("Help Window", "Use mouse/keyboard to interact with plot."\
			+" Press 'h' in plot window to display keyboard shortcuts in terminal.")
		self.SpecPlot.start_event_loop(0)

#This iis the Main program. Pretty simple eh?
if __name__=="__main__":
        app=linefinder_tk(None)
        app.title('LineFinder')#Name of application
        print "Starting TB's Line Finder GUI"
        app.mainloop()
