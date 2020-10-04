"""
Julia workshop 2020
data from
filename, filename2 = "con01_190424_163001.rhs", "con01_RRFD3_s1s1p1p1_M1_000ms_030um_Ch126_f01.jld2"


v1.1 changes how the data is allocated from vector{vector} to a simple matrix
v3: uses pyplot rather than plotly which isn't working (link :x causes errors)!

"""

# directory= "A:/Documents/OneDrive - UNSW/2-Teaching/Julia workshop/"
# cd(directory)

cd(dirname(@__FILE__)) # use this to change directory to location of this script



################   Generate and save data  ##############################
# # run IntanNEO_data_vis_A1x32_scratch and extract data from con01_190424_163001.rhs, then the following
# directory="A:/Documents/OneDrive - UNSW/2-Teaching/Julia workshop/"
# outputdirectory = "A:/Desktop"
# pretrial_time = 0.3 # pre trial time
# trial_t = 1.2 # total trial time in seconds (NB: ensure makes trial_data_time same length as trial_data)
# f_name1=filename2[1:end-5]
# process_Chs="126"  # place your desired channel in quotes
# f_name=string(f_name1[1:38],process_Chs, f_name1[42:end-4])
#
# using FileIO # this is a package that enables you to save the workspace into
# save("$directory$f_name.jld2",
# "fs",fs,
# "pretrial_time",pretrial_time,
# "trial_t",trial_t,
# "laser", laser,
# "f_name", f_name,
# "stimType", stimType,
# "rawSIG", SIG[parse(Int,process_Chs),:])


############################# load data into workspace #########################

######## load data into workspace ########
using Glob , FileIO
data_package=glob("neuro_data.jld2")
rawSIG = load(data_package[1], "rawSIG")
pretrial_time = load(data_package[1], "pretrial_time")
laser = hcat(load(data_package[1], "laser"))
# laser = hcat(load(data_package[1], "laser"))

trial_t1 = load(data_package[1], "trial_t")
fs = load(data_package[1], "fs")

vect_index=collect(1:size(rawSIG,1))  # generate a vector from 1:size of the rawSIG vector

###########

varinfo()  # lets look at what we have in the workspace (variables in Julia's working memory)

#########


size(rawSIG) # this will show you the size of the variable rawSIG

#################
#### plot data for quick inspection

# we can use different plottin libraries:

using PyPlot; pygui()
figure()
# PyPlot.plot(vect_index,rawSIG[1,:],  label = "rawSIG")
# PyPlot.plot(vect_index,5000*laser[1,:], label = "laser displacement")
PyPlot.plot(vect_index,rawSIG,  label = "rawSIG")
PyPlot.plot(vect_index,1000*laser, label = "laser displacement")


using Plots; plotly()
# Plots.plot(vect_index,rawSIG[1,:], hover=vect_index, label = "rawSIG") #plot rawSIG (first row, all columns), display x variable on hover
# Plots.plot!(vect_index,5000*laser[1,:], label = "laser displacement", hover=5000*laser[1,:])

Plots.plot(vect_index,rawSIG, hover=vect_index, label = "rawSIG") #plot rawSIG (first row, all columns), display x variable on hover
Plots.plot!(vect_index,1000*laser, label = "laser displacement", hover=1000*laser[1,:])


# once plotted, you can interact with the plot using the tools above the plot which become visible as you hover your mouse over the plot

##################################

#### zoom in on data region of interest
figure()
# PyPlot.plot(vect_index[209000:212000],rawSIG[1,209000:212000] ,  label = "rawSIG") #plot rawSIG: first row, columns (samples) between 209000 and 212000
# PyPlot.plot(vect_index[209000:212000],5000*laser[1,209000:212000],  label = "laser displacement")

PyPlot.plot(vect_index[209000:212000],rawSIG[209000:212000] ,  label = "rawSIG") #plot rawSIG: first row, columns (samples) between 209000 and 212000
PyPlot.plot(vect_index[209000:212000],1000*laser[209000:212000],  label = "laser displacement")

# plot(vect_index[209000:212000],rawSIG[1,209000:212000] , hover=vect_index, label = "rawSIG") #plot rawSIG: first row, columns (samples) between 209000 and 212000
# plot!(vect_index[209000:212000],5000*laser[1,209000:212000],  label = "laser displacement")

###################################

##### filter the data
using DSP  # this is calling the DSP package which has tools for digital signal processing
bpass = 300 # define Bandpass filter at 330 Hz
bstop = 5000 # define Bandstop filter at 3300 Hz
responsetype = Bandpass(bpass, bstop, fs = fs) # standard filtering use 300-5000
designmethod = Butterworth(6)
# HFdata = filtfilt(digitalfilter(responsetype, designmethod), rawSIG')
HFdata = filtfilt(digitalfilter(responsetype, designmethod), rawSIG)
negHFdata = -HFdata # convention is to invert signal of extracellular recordings (upward deflection represents action potentials)
# HFdata=HFdata' # transpose HFdata matrix (or vector in this case) to be in the same format as rawSIG and laser
varinfo()

#############################################



# view signals:
figure()
# plot(vect_index[209000:212000],rawSIG[1,209000:212000] ,  label = "rawSIG") #plot rawSIG: first row, columns (samples) between 209000 and 212000
# plot(vect_index[209000:212000],5000*laser[1,209000:212000],   label = "laser displacement")
# plot(vect_index[209000:212000],negHFdata[1,209000:212000] ,  label = "negHFdata")

plot(vect_index[209000:212000],rawSIG[209000:212000] ,  label = "rawSIG") #plot rawSIG: first row, columns (samples) between 209000 and 212000
plot(vect_index[209000:212000],1000*laser[209000:212000],   label = "laser displacement")
plot(vect_index[209000:212000],negHFdata[209000:212000] ,  label = "negHFdata")

# plot(vect_index[209000:212000],rawSIG[1,209000:212000] , hover=rawSIG[1,209000:212000], label = "rawSIG") #plot rawSIG: first row, columns (samples) between 209000 and 212000
# plot!(vect_index[209000:212000],5000*laser[1,209000:212000],  hover=vect_index, label = "laser displacement")
# plot!(vect_index[209000:212000],negHFdata[1,209000:212000] , hover=negHFdata[1,209000:212000], label = "negHFdata")


###############################


####### lets chop up our signal based on the peak of each mechanical stimulus
# we need to determine where to split each trial from. Lets base this from our laser signal
# lets find the index for each peak of the laser signal:

using PyCall # import the PyCall library, which allows us to use python commands

scisig = pyimport("scipy.signal") # import the scipy signal processing library from python

stim_ind, pkHeights = scisig.find_peaks(1000*laser[:,1], height = 15, distance = 16000) # doesnt like new syntax
# stim_ind, pkHeights = scisig.find_peaks(1000*laser[1,:], height = 15, distance = 16000)
# use the python find peak function to find all the indexes of peaks that are greater than 90
# only allow one peak per 16000 samples


#######################################

# next, chop up signals according to the indicies of peaks,

# we need do define this function called split trials:
"""
split_trials() takes a single vector and splits it into a 2D array according to the a vector
which contains the indexes of where to split the vector. You define how far before and after
the split location for each trial. This function also permits shifting relative to the split point.

%%%%%%% Inputs %%%%%%%

- amp_signal: is an array containing the data from the electrode of interest

- amp_samplerate: is a value containing the samplerate from recording

- stms: is a row vector that indicates the time that each stimulus artifact occurred in the recording

- trial_length: is the length of the trials that the data will be split into (in samples)

- Lag_time: how long after the stimulus artefact each trial will begin (usually '0', but can be
later for chopping up windows for machine learning etc)

- pre_trial_length: the amount of recording before the stimulus artefact to be included (in samples)

%%%%%%% Outputs %%%%%%%

- trial_data: an n*m matrix conatining all the data for each trial
     n is the number of trials and m is the length of the trials in sample
"""
function split_trials(amp_signal, amp_samplerate, stms, trial_length, lag_time, pre_trial_length)
    trial_data = [];
    # creating the matrices containg each trial as recorded by each of seven electrodes
    for k = 1:size(stms,2)
        strts = convert(Int,floor(stms[1,k]+lag_time)-pre_trial_length)
        t_end = convert(Int,strts+(trial_length-1))
        if t_end > length(amp_signal)
            println("The end of the final trial was too close to the end of the\n",
            "recording to acquire background data for the last signal") &&  break
        else
            # trial_data[k,:] = amp_signal[strts:t_end]
            push!(trial_data, amp_signal[strts:t_end])
        end
    end
    return trial_data
end


# after you define the fuction with the above text you can see it in the help menu (? split_trials)
###############################

# split the HFdata and the laser data into 10 trials:

# trial_data = split_trials(negHFdata[1, :], fs, stim_ind', 5000,  0, 1000)
# trial_laser = split_trials(laser[1, :], fs, stim_ind', 5000,  0, 1000)

trial_data = split_trials(negHFdata, fs, stim_ind', 5000,  0, 1000)
trial_laser = split_trials(laser, fs, stim_ind', 5000,  0, 1000)

# #
# trial_data =  reduce(hcat,[a for a in trial_data ]) # converts array of vectors into a single matrix
# #
# trial_laser = reduce(hcat,[a for a in trial_laser ]) # converts array of vectors into a single matrix
#
# # ml=maximum(length.(trial_data))
# # trial_data2=reduce(hcat,[ if length(a)==ml a else hcat(a,Array{Missing}(missing,ml-length(a))) end for a in trial_data ])
#


#############################################

#############  view signals (single axis)
# first lets convert our vect_index to a time series in ms

time_vect=-1000:1:3999 # this is the 5000 vector starting from -1000
time_vect = (time_vect/fs)*1000


#############  view signals (individual axes)
figure()
PyPlot.plot(time_vect, trial_data[1])
PyPlot.plot(time_vect, 1000*trial_laser[1])
PyPlot.legend(["trial data", "displacement"])
PyPlot.xlabel("milliseconds")
PyPlot.ylabel("trail data (uV) | displacement (um)")

No_trials=size(trial_data,1) # number of trials
fig1, ax1 = subplots(No_trials, 1, sharex="all", sharey="all") # define subplots and which axes are shared
for i = 1:No_trials # for each trail plot the data and the displacement
    ax1[i,1].plot(time_vect, trial_data[i])
    ax1[i,1].plot(time_vect, 1000*trial_laser[i])
    # ax1[i,1].plot((spike_ind[i]/30) .-300, pkH[i],"r.", markersize = 6, color = "red")
end

fig1.suptitle("Signals of $No_trials trials") # add title to subplots
fig1.legend(["trial data", "displacement (microns)"]) # add legend to subplot

############################################

# Plots.plot(time_vect, [trial_data, 5000*trial_laser], layout = (10,1))
# Plots.plot(time_vect, [trial_data, 5000*trial_laser], layout = (10,1), link = :x) #  not working
############################################

###### calculate mean of signals
using Statistics
HFdata_mean=mean(trial_data)


figure()
plot(time_vect,trial_data[1], label = "trial 1")
plot(time_vect,trial_data[5], label = "trial 5")
plot(time_vect,trial_data[10], label = "trial 10")
plot(time_vect,HFdata_mean, color="black", linewidth=2)
legend(["trial 1", "trial 5", "trial 10", "mean of 10 trials"])
# we can see this doesnt capture all the signals well because of the jitter (need to zoom in)

###########################################

# we can see this doesnt capture all the signals well because of the jitter
# lets look for peaks in the first signal

x = trial_data[1] # we will define x and the first row of the matrix trial_data
figure()
plot(time_vect,x) # note I dont need PyPlot.plot, because pyplt was called first so that is the default plotting library
legend(["first trial"])
################################################
# calcualte mean +/- 3 standard deviations to determine what peaks are background noise
# lets calcualte this from the first 1000 samples:
bkgnd_1= 3 * std(trial_data[1][1:1000])

################################################
# now lets find all peaks greater than bkgnd_1, and plot on graph:

ind_x, pks_x = scisig.find_peaks(x, height = bkgnd_1) # use "find_peaks" function from python scisig library to find all peaks in x greater than the height bkgnd_1
plot(time_vect,bkgnd_1*ones(5000,1), label = "3x SD above background noise (bkgnd_1)")
plot((ind_x .-1000)/fs*1000,pks_x["peak_heights"], "r.", markersize = 6, color = "red")
legend(["1st trial", "detection threshold", "peaks above threshold"])
# this last expression gets the index for the peak, subtracts 1000 because we are converting to time and the vector starts at -1000 samples from time zero
# we then divide by the sample frequency (fs) x1000 to convert our x-axis to ms rather than indecies of the vector.


################################################

#lets do this for all trials

inds_trials = [] # need to declare the variable outside a for loop
pks_trials = [] # need to declare the variable outside a for loop
for i = 1: 10
    ind_x, pks_x = scisig.find_peaks(trial_data[i], height = 3 * std(trial_data[i][1:1000]))
    push!(inds_trials, ind_x')
    push!(pks_trials, pks_x["peak_heights"]')
end
################################################

# test it works on a randon example:

figure()
plot(time_vect,trial_data[7])
plot((inds_trials[7] .-1000)[1,:]/fs*1000,pks_trials[7], "r.", markersize = 6, color = "red")
# looks ok; try substituting another number instead of 7...
# can anyone tell me why we need the [1,:]?

################################################
# lets plot all the signals on separate subplots:

fig2, ax2 = subplots(No_trials, 1, sharex="all", sharey="all") # define subplots and which axes are shared
for i = 1:No_trials # for each trail plot the data and the displacement
    ax2[i,1].plot(time_vect, trial_data[i])
    ax2[i,1].plot(time_vect, 1000*trial_laser[i])
    ax2[i,1].plot((inds_trials[i] .-1000)[1,:]/fs*1000, pks_trials[i],"r.", markersize = 6, color = "red")
end
fig2.suptitle("Signals of $No_trials trials") # add title to subplots
fig2.legend(["trial data", "displacement (microns)", "detected peaks"]) # add legend to subplot


################################################
# Lets make a raster for all 10 trials and see what it looks like

figure()
for i =1:No_trials
plot((inds_trials[i] .-1000)/fs*1000, i*ones(size(pks_trials[i])),"r.", markersize = 6, color = "red")
end
# note that the raster is in the reverse order to the previous figure; why?
################################################

# Lets add a frequency histogram to count the number of events above threshold:

binwidth=1 #define bin width
pretrial_length = 1000 # each trial has 1000 samples before it
trial_len = 5000 # each trial has 5000 samples in total
bin_edges=collect(round(-pretrial_length/fs*1000):binwidth: round((-pretrial_length+trial_len)/fs*1000))

inds_trials_flat = sort(reduce(hcat,inds_trials)') # concatinate all vectors horixontally, then sort in order

inds_trials_flat = inds_trials_flat .-1000 # subtract 1000 (range = -1000:4000, i.e. 5000 length vector)

hist(inds_trials_flat/fs*1000, bins=bin_edges, alpha = 0.5)

## now add to the raster and link the axes
fig3, ax3 = subplots(2, 1, sharex="all") # start a new plot
for i =1:No_trials
ax3[1,1].plot((inds_trials[i] .-1000)/fs*1000, i*ones(size(pks_trials[i])),"r.", markersize = 6, color = "red")
end
ax3[2,1].hist(inds_trials_flat/fs*1000, bins=bin_edges, alpha = 0.5)
