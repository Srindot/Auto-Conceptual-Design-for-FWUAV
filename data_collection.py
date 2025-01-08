from simulation import Mk4SaveDataToCSV, simulation
import pterasoftware as ps
import numpy as np
import matplotlib.pyplot as plt
from SelfFunctions import extract_forces, extract_second_cycle
import os 
import csv



def Mark4Simulation(FlappingPeriods, Angles_of_Attacks, Air_Speeds, mw_wingspans, mw_root_chords, taper_ratio):
    n = 0
    
    # Function call
    for w in  mw_wingspans:
        for e in mw_root_chords:
            for r in taper_ratio:
                for i in FlappingPeriods:
                    for j in Air_Speeds:
                        for k in Angles_of_Attacks:
                                
                                lift, induced_drag  = simulation(mw_airfoil = "naca8304", fp = i, va = j, aoa = k, 
                                                                                 mw_wingspan = w, mw_root_chord = e, taper_ratio = r)
                                # Flatten the arrays
                                lift = lift.flatten()
                                induced_drag = induced_drag.flatten()                                
                                print("Before Second Cycle Extraction :  ", lift.shape, induced_drag.shape)
                                print("Before Second Cycle Extraction value:  ", lift, induced_drag)

                                # Extract the second cycle from the data
                                lift, induced_drag,  = extract_second_cycle(lift, induced_drag)
                                print("After Second Cycle Extraction :  ",lift.shape, induced_drag.shape)
                                print("After Second Cycle Extraction value:  ", lift, induced_drag )

                                # Save data to CSV with a unique filename for each combination
                                file_name = f"Data/Data_Instance{n}.csv"
                                Mk4SaveDataToCSV("naca8304", i, j, k, w, e, r, lift, induced_drag, file_name)
                                n += 1
 
    print("Data Collection is Over")



FlappingPeriods = [0.4]
Angles_of_Attacks = [5]
Air_Speeds =  [5]
mw_wingspans = [1.5]
mw_root_chords = [0.3]
taper_ratio = [0.66]
Mark4Simulation(FlappingPeriods, Angles_of_Attacks, Air_Speeds, mw_wingspans, mw_root_chords, taper_ratio)