# LineFinder
A Python GUI to identify spectral lines of systems at different redshifts/velocities.

Author: Trystyn Berg

This readme contains a brief outline of the purpose of the code, the input/output requirements and formats, and a quick start guide to using the code.

This GUI was written for astronomers who are identifying many quasar absorption line systems in the same spectrum. However, the code can be used for any other purposes where a given spectrum is shifted by a given velocity and requires identifying which lines might be present at those velocities.

The handy feature of this program is to keep track of all the diferent absorption systems (i.e. different velocities). This is done by saivng all the information to a logfile. Built into the logfile, the user can add velocity limits on the absorption profile, the colour of the line (for display purposes), and a flag to say if the line is saturated, blended, good, etc. This may be useful if you would like to calculate abundances later. The format of this logfile is given below.

PYTHON PACKAGE REQUIREMENTS:

       	TB_LINEFINDER was developed with:
       	
               	-Python Version 2.7
               	
               	-Matplotlib Version 1.4.1
               	
                       	-Need matplotlib.backends.backend_tkagg
                       	
               	-Tkinter Revision: 81008
               	
               	-Numpy
               	
       	Other verisons of above packages may work as well.
       	
       	
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

OUTPUT LOGFILE FORMAT:

       	-The output is a semicolon (;) delimited
  
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

                     The recommended FLAG notation is a binary format with the following option
  
                            0 - No good/skip
       
                            1 - OK
       	
                            2 - Blend
       	
                            4 - Upper Limit/Non-detection
       
                            8 - Lower Limit/Saturated

                     If one prefers a different format, the user can use any integer. FLAG is only used for the purposes of saving to the log as a reminder


                      When the output logfile is being read, any line beginning with '#' will be ignored

                     For the purpose of future reference, the first two lines of the logfile contain the input filenames in the format:

                             #!Line List: <INPUT LINELIST ASCII FILENAME>

                            #!Spectrum: <INPUT SPECTRUM ASCII FILENAME>



QUICK USEAGE:

       Once you have created the necessary input files, start-up tb_linefinder. Locate the files you created on the hard drive, and name the output logfile. If you so desire, you can load the inputs/output via command line by calling the program as:
       
       tb_linefinder.py '<INPUT_SPECTRUM_FILENAME>' '<INPUT_LINELIST_FILENAME>' '<OUTPUT_LOGFILE_NAME>'
    
       A tutorial pop-up window is enabled by default. You can disable either in the main menu or completely via the code.
    
       Once the FIND LINES button is pressed, the spectrum will be plotted. Use the mouse and keyboard to indtify lines:
    
             -placing the cursor on a given feature and typing 'l' will assume the feature is Lyman alpha and calcualte the redshift. From there it will display the velocity profiles of all the lines at that redshift. The user then selects which lines are there, and can then add notes (via the terminal window) to the log file.
      
             -Similar to 'l', placing the cursor on a given feature and typing 'a' will allow you to choose the line (via the terminal) and caculate the redshift. Again all the velocity profiles are plotted, and you select which are real.\
      
             -You can also hit 'z', and type in the redshift into the terminal window. Selecting lines will happen the same way as hitting 'l' or 'a'
      
             - You can edit the logfile by pressing 'e'. A drop down menu and buttons allows you to interact. If you click the EDIT LOGFILE button, it will interact with you via the terminal window.
       
             -'q' will quit the linefinder routine and take you back to the main window.
      
             -If you ever forget the keys, hit 'h', and a list of options show up!
