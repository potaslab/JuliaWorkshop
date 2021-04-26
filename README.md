# JuliaWorkshop
Content for getting started with Julia.

This repository contains the following files: 

*	"Julia workshop 2020 slides github.pdf" which explains some basics on data acquisition and how the data has been captured

*	“neuro_data.jld2” which contains the data

*	“Installing Julia and Jupyter 2020” which contains instructions to get you started with Julia using Jupyter notebook

*	“JuliaWorkshop2020.ipynb” which is a Jupyter notebook file. 

*	“JuliaWorkshop2020_v3.jl” which is a Julia script with similar content as the .ipynb file.

Note that there are no instructions in the attached files for the .jl file as yet. For this  you will need to install Atom and Juno.  


# Instructions for Atom: 

*	Once you have installed Julia on  your computer, install Atom for you operating system: https://atom.io/

*	Once installed, go to the settings (Ctrl + ,) -> intall -> search for “uber-juno”

*	Other worthwhile packages to install are: auto-indent, copy-filename, highlight-selected, simple-drag-drop-text


# Installing Julia/VS Code/VS Code Julia extension

* Install Julia for your platform: https://julialang.org/downloads/

* Install VS Code for your platform: https://code.visualstudio.com/download. At the end of this step you should be able to start VS Code.

* Install the Julia VS Code extension:

*   Start VS Code.

*   Inside VS Code, go to the extensions view either by executing the *View: Show Extensions* command (click View->Command Palette...) or by clicking on the extension icon on the left side of the VS Code window.

*   In the extensions view, simply search for the term *julia* in the marketplace search box, then select the extension named *Julia* and click the install button. You might have to restart VS Code after this step.



# Instructions to update Julia:

To update Julia, you must check that Julia-client is correctly configured:


*	Settings -> packages -> Julia-client -> settings -> Julia Path: 

C:\Users\userprofile\AppData\Local\Programs\Julia-1.6.1\bin\julia.exe

Ensure that it correctly points to the Julia.exe as shown above. 

# Exercise: 
Find the number of events evoked in the first 25 ms for the 8th trail under two different bandpass filtering conditions, using a detection threshold of 3 standard deviations above background activity: 

* 300 - 5000 Hz

* 1000 - 2000 Hz

Plot these two filtered signals on a single figure, showing the detection threshold and identified peaks above threshold. Ensure your figure is fully labelled (x and y axes, figure legend and title etc).
