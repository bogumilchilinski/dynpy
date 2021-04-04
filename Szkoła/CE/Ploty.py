import matplotlib.pyplot as plt 
import numpy as np



    


def plot(x,y,marker = None):

    fig = plt.figure(figsize = (15,10) , dpi=80)
    ax = plt.subplot()
    ax.grid(True)
    ax.set_axisbelow(True)        

    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)


    return ax.plot(x, y, marker = marker)





def plot_double(x,y1,y2):


    color1 = 'black'
    color2 = 'red'

    fig = plt.figure(figsize = (15,10) , dpi=80)
    ax = plt.subplot()
    ax2 = ax.twinx()
    ax.grid(True)
    ax.set_axisbelow(True)

    #ax.set_xlabel('Engine speed [RPM]', fontsize = 20)
    ax.set_ylabel('Torque [Nm]',fontsize = 20, color='black')
    ax2.set_ylabel('Power [kW]',fontsize = 20 , color='red')
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='y', labelsize=20)
    ax.set_ylim(0,280)
    ax2.set_ylim(0, 150)

    

    return ax.plot(x, y1, color=color1) , ax2.plot(x, y2, color=color2) , ax.scatter(x,y1,color=color1), ax2.scatter(x,y2, color=color2)

        








        






