#!/usr/bin/env python3

import numpy as np
import logging
import os.path
import time
import math

# logger setup
logger = logging.getLogger(__name__)

##### USER DEFINED GENERAL SETTINGS #####

#set new name for each experiment, otherwise files will be overwritten
EXP_NAME = 'cas9pace_20210217_darwin_expt'
EVOLVER_IP = '192.168.1.9'
EVOLVER_PORT = 8081

##### Identify pump calibration files, define initial values for temperature, stirring, volume, power settings

TEMP_INITIAL = [37] * 16 #degrees C, makes 16-value list
#Alternatively enter 16-value list to set different values
#TEMP_INITIAL = [30,30,30,30,32,32,32,32,34,34,34,34,36,36,36,36]

STIR_INITIAL = [8] * 16 #try 8,10,12 etc; makes 16-value list
#Alternatively enter 16-value list to set different values
#STIR_INITIAL = [7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10]

LAGOON_VOLUME = 10 # ml
VOLUME =  30 #mL, determined by vial cap straw length
PUMP_CAL_FILE = 'pump_cal.txt' #tab delimited, mL/s with 16 influx pumps on first row, etc.
OPERATION_MODE = 'chemostat' #use to choose between 'turbidostat' and 'chemostat' functions


# IPP Calibrations. Function of form rate = c*frequency^b
c = [5.016]
b = [.5485]

efflux_addrs = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 ,30, 31]

##### END OF USER DEFINED GENERAL SETTINGS #####

def hz_to_rate(hz, c, b):
    return c * math.pow(hz, b)

def rate_to_hz(rate, c, b):
    """ rate should be in ml/h """
    return math.pow(((rate) / c), 1.0/b)

def ipp_calculations(elapsed_time, eVOLVER):

    # arabinose stock concentration is at 250 mM
    # Solenoid addresses for each ipp.
    # Each pump requires 3 addresses. These vars capture the 1st of the three (sequential)
    v2v_addr = 32
    ipp_min_waiting_time = 4

    # Sets the minimum amount of time that the experiment must run
    # before the IPP selection scheme will start
    bolus_amount = .4 # ml
    bolus_rate = 10
    bolus_time = bolus_amount / hz_to_rate(10, c[0], b[0])

    turnover_time = 1 # hours
    init_rate = .08
    start_rate = 0.04 # V/h

    # Start ipp protocol
    if (elapsed_time > ipp_min_waiting_time):
        if (elapsed_time < ipp_min_waiting_time + bolus_time):
            # Start bolus      
            print("ipp bolus")      
            rate = bolus_rate
            vial = 'all'
        else:
            rate = rate_to_hz(start_rate * LAGOON_VOLUME, c[0], b[0])
            vial = 'all'

        print("running ipp cmd. addr: {0}, vial: {1}, rate: {2}".format(v2v_addr, vial, round(rate,3)))
        eVOLVER.ipp_command(v2v_addr, vial, round(rate,3))    

def turbidostat(eVOLVER, input_data, vials, elapsed_time, run_efflux):
    OD_data = input_data['transformed']['od']

    ##### USER DEFINED VARIABLES #####

    turbidostat_vials = vials #vials is all 16, can set to different range (ex. [0,1,2,3]) to only trigger tstat on those vials
    stop_after_n_curves = np.inf #set to np.inf to never stop, or integer value to stop diluting after certain number of growth curves

    #Alternatively, use 16 value list to set different thresholds, use 9999 for vials not being used
    lower_thresh = [0.5, 0.5, 0.5, 0.5, 9999, 9999, 9999, 9999, 0.5, 0.4, 0.5, 0.5, 9999, 9999, 9999, 9999]
    upper_thresh = [0.6, 0.6, 0.6, 0.6, 9999, 9999, 9999, 9999, 0.6, 0.5, 0.6, 0.6, 9999, 9999, 9999, 9999]

    ##### END OF USER DEFINED VARIABLES #####


    ##### Turbidostat Settings #####
    #Tunable settings for overflow protection, pump scheduling etc. Unlikely to change between expts

    time_out = 5 #(sec) additional amount of time to run efflux pump
    pump_wait = 3 # (min) minimum amount of time to wait between pump events

    ##### End of Turbidostat Settings #####

    save_path = os.path.dirname(os.path.realpath(__file__)) #save path
    flow_rate = eVOLVER.get_flow_rate() #read from calibration file


    ##### Turbidostat Control Code Below #####

    #ipp_calculations(elapsed_time, eVOLVER)

    # fluidic message: initialized so that no change is sent
    MESSAGE = ['--'] * 48
    if run_efflux:
        for addr in efflux_addrs:
            MESSAGE[addr] = 6
    for x in turbidostat_vials: #main loop through each vial

        # Update turbidostat configuration files for each vial
        # initialize OD and find OD path

        file_name =  "vial{0}_ODset.txt".format(x)
        ODset_path = os.path.join(save_path, EXP_NAME, 'ODset', file_name)
        data = np.genfromtxt(ODset_path, delimiter=',')
        ODset = data[len(data)-1][1]
        ODsettime = data[len(data)-1][0]
        num_curves=len(data)/2;

        file_name =  "vial{0}_OD.txt".format(x)
        OD_path = os.path.join(save_path, EXP_NAME, 'OD', file_name)
        data = np.genfromtxt(OD_path, delimiter=',')
        average_OD = 0

        # Determine whether turbidostat dilutions are needed
        enough_ODdata = (len(data) > 7) #logical, checks to see if enough data points (couple minutes) for sliding window
        collecting_more_curves = (num_curves <= (stop_after_n_curves + 2)) #logical, checks to see if enough growth curves have happened

        if enough_ODdata:
            # Take median to avoid outlier
            od_values_from_file = []
            for n in range(1,7):
                od_values_from_file.append(data[len(data)-n][1])
            average_OD = float(np.median(od_values_from_file))

            #if recently exceeded upper threshold, note end of growth curve in ODset, allow dilutions to occur and growthrate to be measured
            if (average_OD > upper_thresh[x]) and (ODset != lower_thresh[x]):
                text_file = open(ODset_path, "a+")
                text_file.write("{0},{1}\n".format(elapsed_time,
                                                   lower_thresh[x]))
                text_file.close()
                ODset = lower_thresh[x]
                # calculate growth rate
                eVOLVER.calc_growth_rate(x, ODsettime, elapsed_time)

            #if have approx. reached lower threshold, note start of growth curve in ODset
            if (average_OD < (lower_thresh[x] + (upper_thresh[x] - lower_thresh[x]) / 3)) and (ODset != upper_thresh[x]):
                text_file = open(ODset_path, "a+")
                text_file.write("{0},{1}\n".format(elapsed_time, upper_thresh[x]))
                text_file.close()
                ODset = upper_thresh[x]

            #if need to dilute to lower threshold, then calculate amount of time to pump
            if average_OD > ODset and collecting_more_curves:

                time_in = - (np.log(lower_thresh[x]/average_OD)*VOLUME)/flow_rate[x]

                if time_in > 20:
                    time_in = 20

                time_in = round(time_in, 2)

                file_name =  "vial{0}_pump_log.txt".format(x)
                file_path = os.path.join(save_path, EXP_NAME,
                                         'pump_log', file_name)
                data = np.genfromtxt(file_path, delimiter=',')
                last_pump = data[len(data)-1][0]
                if ((elapsed_time - last_pump)*60) >= pump_wait: # if sufficient time since last pump, send command to Arduino
                    logger.info('turbidostat dilution for vial %d' % x)
                    # influx pump
                    MESSAGE[x] = str(time_in)
                    # efflux pump
                    MESSAGE[x + 16] = str(time_in + time_out)

                    file_name =  "vial{0}_pump_log.txt".format(x)
                    file_path = os.path.join(save_path, EXP_NAME, 'pump_log', file_name)

                    text_file = open(file_path, "a+")
                    text_file.write("{0},{1}\n".format(elapsed_time, time_in))
                    text_file.close()
        else:
            logger.debug('not enough OD measurements for vial %d' % x)

    # send fluidic command only if we are actually turning on any of the pumps
    if MESSAGE != ['--'] * 48:
        eVOLVER.fluid_command(MESSAGE)

        # your_FB_function_here() #good spot to call feedback functions for dynamic temperature, stirring, etc for ind. vials
    # your_function_here() #good spot to call non-feedback functions for dynamic temperature, stirring, etc.

    # end of turbidostat() fxn

def chemostat(eVOLVER, input_data, vials, elapsed_time, run_efflux):
    OD_data = input_data['transformed']['od']

    ##### USER DEFINED VARIABLES #####
    start_OD = 0 # ~OD600, set to 0 to start chemostate dilutions at any positive OD
    start_time = 4 #hrs, set 0 to start immediately
    # Note that script uses AND logic, so both start time and start OD must be surpassed

    chemostat_vials = vials #vials is all 16, can set to different range (ex. [0,1,2,3]) to only trigger tstat on those vials

    #rate_config = [1] * 16 #to set all vials to the same value, creates 16-value list
    #UNITS of 1/hr, NOT mL/hr, rate = flowrate/volume, so dilution rate ~ growth rate, set to 0 for unused vials

    #Alternatively, use 16 value list to set different rates, use 0 for vials not being used

    # SLOW PUMP ARRAY FOR PACE: First 8 (0-7) pumps are for media/chemostat, last 8 (7-15) are for vial to vial.
    #rate_config = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    rate_config = [1,1,1,1,.25,1,1.25,.5,1,1,1,1,.25,1.25,.75,.75]

    ##### END OF USER DEFINED VARIABLES #####


    ##### Chemostat Settings #####

    #Tunable settings for bolus, etc. Unlikely to change between expts
    bolus = 0.1 #mL, slow rates can do much lower bolus sizes

    ##### End of Chemostat Settings #####

    save_path = os.path.dirname(os.path.realpath(__file__)) #save path
    flow_rate = eVOLVER.get_flow_rate() #read from calibration file
    period_config = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] #initialize array
    bolus_in_s = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] #initialize array


    ##### Chemostat Control Code Below #####

    for x in chemostat_vials: #main loop through each vial

        # Update chemostat configuration files for each vial

        #initialize OD and find OD path
        file_name =  "vial{0}_OD.txt".format(x)
        OD_path = os.path.join(save_path, EXP_NAME, 'OD', file_name)
        data = np.genfromtxt(OD_path, delimiter=',')
        average_OD = 0
        enough_ODdata = (len(data) > 7) #logical, checks to see if enough data points (couple minutes) for sliding window

        if enough_ODdata: #waits for seven OD measurements (couple minutes) for sliding window

            #calculate median OD
            od_values_from_file = []
            for n in range(1, 7):
                od_values_from_file.append(data[len(data)-n][1])
            average_OD = float(np.median(od_values_from_file))

            # set chemostat config path and pull current state from file
            file_name =  "vial{0}_chemo_config.txt".format(x)
            chemoconfig_path = os.path.join(save_path, EXP_NAME,
                                            'chemo_config', file_name)
            chemo_config = np.genfromtxt(chemoconfig_path, delimiter=',')
            last_chemoset = chemo_config[len(chemo_config)-1][0] #should t=0 initially, changes each time a new command is written to file
            last_chemophase = chemo_config[len(chemo_config)-1][1] #should be zero initially, changes each time a new command is written to file
            last_chemorate = chemo_config[len(chemo_config)-1][2] #should be 0 initially, then period in seconds after new commands are sent

            # once start time has passed and culture hits start OD, if no command has been written, write new chemostat command to file
            #if ((elapsed_time > start_time) and average_OD > start_OD):
            if ((elapsed_time > start_time)):

                #calculate time needed to pump bolus for each pump
                bolus_in_s[x] = bolus/flow_rate[x]

                # calculate the period (i.e. frequency of dilution events) based on user specified growth rate and bolus size
                if rate_config[x] > 0:
                    c_vials = [0,1,2,3,8,9,10,11]
                    if x in c_vials:
                        period_config[x] = (3600*bolus)/((rate_config[x])*VOLUME) #scale dilution rate by bolus size and volume
                    else:
                        period_config[x] = (3600*bolus)/((rate_config[x])*LAGOON_VOLUME) #scale dilution rate by bolus size and volume
                else: # if no dilutions needed, then just loops with no dilutions
                    period_config[x] = 0

                if  (last_chemorate != period_config[x]):
                    print('Chemostat updated in vial {0}'.format(x))
                    logger.info('chemostat initiated for vial %d, period %.2f'
                                % (x, period_config[x]))
                    # writes command to chemo_config file, for storage
                    text_file = open(chemoconfig_path, "a+")
                    text_file.write("{0},{1},{2}\n".format(elapsed_time,
                                                           (last_chemophase+1),
                                                           period_config[x])) #note that this changes chemophase
                    text_file.close()
        else:
            logger.debug('not enough OD measurements for vial %d' % x)

        # your_FB_function_here() #good spot to call feedback functions for dynamic temperature, stirring, etc for ind. vials
    # your_function_here() #good spot to call non-feedback functions for dynamic temperature, stirring, etc.

    eVOLVER.update_chemo(input_data, chemostat_vials, bolus_in_s, period_config) #compares computed chemostat config to the remote one
    time.sleep(1)
    ipp_calculations(elapsed_time, eVOLVER)
    time.sleep(1)
    MESSAGE = ['--'] * 48
    if run_efflux:
        for addr in efflux_addrs:
            MESSAGE[addr] = 15
        eVOLVER.fluid_command(MESSAGE)
    # end of chemostat() fxn

    eVOLVER.update_chemo(input_data, chemostat_vials, bolus_in_s, period_config)

# def your_function_here(): # good spot to define modular functions for dynamics or feedback

if __name__ == '__main__':
    print('Please run eVOLVER.py instead')
