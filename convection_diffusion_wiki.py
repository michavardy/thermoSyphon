import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
import math
from tqdm import tqdm
from iapws import IAPWS97

class Constant(Enum):

    #### material Constants
    SPECIFIC_HEAT = 4184 ## Joule / Kilogram * Kalvin
    MASS_DENSITY = 997 #997 ## Kilogram / meters ^3
    THERMAL_CONDUCTIVITY = 0.598 ## Watt / meter Kalvin
    DIAMETER_CONTAINER = 0.3 ## Meters 0.3, 0.07
    HEIGHT_CONTAINER = 0.3 ## Meters 0.3, 0.07
    AREA_CONTAINER_SECTION = math.pi * (DIAMETER_CONTAINER / 2)**2 ## Meters ^2
    HOT_TEMPERATURE = 70 + 273.15 ## Kalvin
    INITIAL_TEMPERATURE = 25 + 273.15  ## Kalvin
    HEAT_SOURCE = 6000 ## Watt 5000
    MASS_RATE = HEAT_SOURCE / (SPECIFIC_HEAT * (HOT_TEMPERATURE - INITIAL_TEMPERATURE))  ## Kilograms / Second
    VELOCITY =  MASS_RATE / (MASS_DENSITY * AREA_CONTAINER_SECTION)  ## Meters / Second
    DIFFUSIVITY =  THERMAL_CONDUCTIVITY / (SPECIFIC_HEAT * MASS_DENSITY) ## Meters ^2 / Second


    #### Numerical Constants
    DURATION = 80 * 60 ## Second
    #TOTAL_NUMBER_OF_SPACE_STEPS = 1000 # 1000
    #TOTAL_NUMBER_OF_TIME_STEPS = 10000 # 10000
    #SPACE_STEP = HEIGHT_CONTAINER / TOTAL_NUMBER_OF_SPACE_STEPS  
    #TIME_STEP = DURATION / TOTAL_NUMBER_OF_TIME_STEPS ## Second
    SPACE_STEP = 0.003 # meters 0.0003
    TIME_STEP = 0.06 # Seconds 0.06
    TOTAL_NUMBER_OF_SPACE_STEPS = int(HEIGHT_CONTAINER / SPACE_STEP)
    TOTAL_NUMBER_OF_TIME_STEPS = int(DURATION / TIME_STEP)



    

    #### INITIALIZE 
    INTIAL_TEMPERATURE_DISTRIBUTION_LIST_OVER_SPACE_DOMAIN = [HOT_TEMPERATURE] + [INITIAL_TEMPERATURE] * (TOTAL_NUMBER_OF_SPACE_STEPS - 1) ### list of floats with initial temp dist
class Output:
    name = 'simulation_results'
    path = f'{name}.xlsx'
    df = pd.DataFrame()
class SolveCoeff:
    def __init__(self):
        pass
    def const():
        """
            not so accurate
        """
        DIFF = ( Constant.DIFFUSIVITY.value * Constant.TIME_STEP.value) / (Constant.SPACE_STEP.value**2) ## DIFFUSION TERM
        ADV = (Constant.VELOCITY.value * Constant.TIME_STEP.value) / (2 * Constant.SPACE_STEP.value) ## ADVECTION (flow) TERM
        A = 1 - 2 * DIFF
        B = DIFF + ADV
        C = DIFF - ADV
        return (A,B,C)
    def nonLin(temp):
        """
            this one takes forever, never use it
        """
        MASS_DENSITY = IAPWS97(T=temp, P=0.101325).rho
        DIFFUSIVITY = Constant.THERMAL_CONDUCTIVITY.value / (Constant.SPECIFIC_HEAT.value * MASS_DENSITY) ## Meters ^2 / Second
        DIFF = (DIFFUSIVITY * Constant.TIME_STEP.value) / (Constant.SPACE_STEP.value**2) ## DIFFUSION TERM
        ADV = (Constant.VELOCITY.value * Constant.TIME_STEP.value) / (2 * Constant.SPACE_STEP.value) ## ADVECTION (flow) TERM
        A = 1 - 2 * DIFF
        B = DIFF + ADV
        C = DIFF - ADV
        return (A,B,C)
    def weighted():
        """
            in testing
        """
        DIFF = 4 * ( Constant.DIFFUSIVITY.value * Constant.TIME_STEP.value) / (Constant.SPACE_STEP.value**2) ## DIFFUSION TERM
        ADV = 4 * (Constant.VELOCITY.value * Constant.TIME_STEP.value) / (2 * Constant.SPACE_STEP.value) ## ADVECTION (flow) TERM
        A = 1 - 2 * DIFF
        B = DIFF + ADV
        C = DIFF -  ADV
        return (A,B,C)
    def cgpt():
        b =  0.2*Constant.VELOCITY.value #1.5
        D =  1 * Constant.DIFFUSIVITY.value #3
        A = (D*Constant.TIME_STEP.value/(Constant.SPACE_STEP.value**2)) 
        B = (b*Constant.TIME_STEP.value/Constant.SPACE_STEP.value)
        return A,B
def check_stability():
    s1 =  Constant.SPACE_STEP.value < ((2 * Constant.DIFFUSIVITY.value) / Constant.VELOCITY.value)
    s2 = Constant.TIME_STEP.value < ((Constant.SPACE_STEP.value **2) / (2* Constant.DIFFUSIVITY.value))
    print(f"Numerical Constant Evaluation: Time Step: {Constant.TIME_STEP.value}, Space Step: {Constant.SPACE_STEP.value}, velocity: {Constant.VELOCITY.value}")
    print(f"Diffusivity: {Constant.DIFFUSIVITY.value}")
    print(f"Stability Evaluation: Space Step Stability: {s1}, Time Step Stability: {s2}")

def solve_convection_diffusion():
    # point temp_dist to initial temperature distribution over time
    all_temp_dist = [Constant.INTIAL_TEMPERATURE_DISTRIBUTION_LIST_OVER_SPACE_DOMAIN.value.copy()]

    # loop over total number of time steps starting at 2nd time step
    for time_step in tqdm(range(1, Constant.TOTAL_NUMBER_OF_TIME_STEPS.value), desc="loop over total number of time steps"):
        
        # we take the temperature distribution from the previous time step
        temp_dist_prev = all_temp_dist[time_step - 1]

        # we take the temperature distribution from the current time step
        temp_dist_cur = Constant.INTIAL_TEMPERATURE_DISTRIBUTION_LIST_OVER_SPACE_DOMAIN.value.copy()

        # loop over total number of space steps starting at 2nd space step
        for space_step in range(1, Constant.TOTAL_NUMBER_OF_SPACE_STEPS.value - 1):
            # calculate coefficients
            A,B,C = SolveCoeff.weighted()
            # solve the current time step using the previous time step
            temp_dist_cur[space_step] = temp_dist_prev[space_step] * A + temp_dist_prev[space_step - 1] * B + temp_dist_prev[space_step + 1] * C
        # append the current temperature distribuition to all_temp_dist
        all_temp_dist.append(temp_dist_cur)
    # return all temp dist
    return all_temp_dist

def solve_convection_diffusion_cgpt():
    # point temp_dist to initial temperature distribution over time
    all_temp_dist = [Constant.INTIAL_TEMPERATURE_DISTRIBUTION_LIST_OVER_SPACE_DOMAIN.value.copy()]

    # loop over total number of time steps starting at 2nd time step
    for time_step in tqdm(range(1, Constant.TOTAL_NUMBER_OF_TIME_STEPS.value), desc="loop over total number of time steps"):
        
        # we take the temperature distribution from the previous time step
        temp_dist_prev = all_temp_dist[time_step - 1]

        # we take the temperature distribution from the current time step
        temp_dist_cur = Constant.INTIAL_TEMPERATURE_DISTRIBUTION_LIST_OVER_SPACE_DOMAIN.value.copy()

        # loop over total number of space steps starting at 2nd space step
        for space_step in range(1, Constant.TOTAL_NUMBER_OF_SPACE_STEPS.value - 1):
            # calculate coefficients
            A,B = SolveCoeff.cgpt()
            # solve the current time step using the previous time step
            temp_dist_cur[space_step] = temp_dist_prev[space_step] + A*(temp_dist_prev[space_step + 1] - 2*temp_dist_prev[space_step] + temp_dist_prev[space_step - 1]) - B*(temp_dist_prev[space_step] - temp_dist_prev[space_step - 1]) 
        # append the current temperature distribuition to all_temp_dist
        all_temp_dist.append(temp_dist_cur)
    # return all temp dist
    return all_temp_dist

def isolate_point_over_time(temp_dist, depth_percent):
    index = int(Constant.TOTAL_NUMBER_OF_SPACE_STEPS.value * depth_percent / 100)
    temp_point = [temp[index] - 273.15 for temp in temp_dist]
    return temp_point

def format_temperature_array(array):
    duration_minutes = Constant.DURATION.value / 60
    minutes_array = np.linspace(0, duration_minutes, Constant.TOTAL_NUMBER_OF_TIME_STEPS.value)[::100]
    array = np.asarray(array)[::100]
    return array, minutes_array

def writeToOutput(time,temp,col_name):
    if output.df.empty:
        output.df.time_col = pd.Series(time)
    else:
        output.df[col_name] = pd.Series(temp) 

def plot_points_at_depth(temp_dist, depth_list):
    for depth in depth_list:
        temp_dist_at_point = isolate_point_over_time(temp_dist, depth)
        temp_dist_at_point, time_array = format_temperature_array(temp_dist_at_point)
        writeToOutput(time_array, temp_dist_at_point,f'depth_{depth}')
        plt.plot(time_array, temp_dist_at_point, label=f"depth: {depth}%")
    # format plot
    plt.xlabel('time [minutes]')
    plt.ylabel('Temperature [C]')
    plt.title('Convection Diffusion Water Heating Thermosyphone')
    plt.xlim(-0.5, 80)
    plt.legend()
    plt.show()

def main():
    global output
    output = Output()
    check_stability()
    temp_dist = solve_convection_diffusion_cgpt()
    plot_points_at_depth(temp_dist, [1,50,90])
    output.df.to_excel(output.path, index=False)


    #print(temp_dist[-1])
    
if __name__ == "__main__":
    main()





