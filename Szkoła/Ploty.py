import matplotlib.pyplot as plt 
import numpy as np



    


def plot(x,y,marker = None, label = None):
    fig = plt.figure(figsize = (10,8) , dpi=80)
    ax = plt.subplot()
    ax.grid(True)
    ax.set_axisbelow(True)        

    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

    return ax.plot(x, y, marker = marker, label = label)





def plot_double(x,y1,y2):
    color1 = 'black'
    color2 = 'red'

    fig = plt.figure(figsize = (10,8) , dpi=80)
    ax = plt.subplot()
    ax2 = ax.twinx()
    ax.grid(True)
    ax.set_axisbelow(True)

    ax.set_xlabel('x, Max. speed/base speed', fontsize = 20)
    ax.set_ylabel('Tractive power [kW]',fontsize = 20, color='black')
    ax2.set_ylabel('$ dP_t/dV_b  ~ [kW/(m/s^2)] $',fontsize = 20 , color='red')
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='y', labelsize=20)
    #ax2.set_ylim(-220,220)
    #ax.set_ylim(-8, 8)
    #ax.set_xlim(0, 140000) 
    return ax.plot(x, y1, color=color1) , ax2.plot(x, y2, color=color2) , ax.scatter(x,y1,color=color1), ax2.scatter(x,y2, color=color2)

        








        






