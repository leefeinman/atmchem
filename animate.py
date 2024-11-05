#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 17:21:24 2024

@author: leefeinman
"""
#%% Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

#%% Function
def animate(line, fig, ax, filename = None, length_secs = 5, set_margins = False, margin_pcnt = 0.05):
    '''
    Animate a subplot.


    Parameters
    ----------
    line : A 'plot' object.
        The line to be animated..
    fig : A fig object.
        DESCRIPTION.
    ax : An axis object.
        DESCRIPTION.
    filename : string, optional
        If you want to save the animation, provide the filename. .GIF is a common filetype. The default is None.
    length_secs : int, optional
        Duration of the animation in seconds. The default is 5.
    set_margins : boolean, optional
        If the xy data being animated should set the xy axis limits, then make this value True. The default is False. 
    margin_pcnt : float, optional
        The scaling margin amount on the axes. Only employed when set_margins is True. The default is 0.05.

    Returns
    -------
    animation.FuncAnimation object.
        

    '''
    line = line[0]
    
    x_vals = line.get_xdata()
    y_vals = line.get_ydata()
    
    # axis limits
    if set_margins:
        ax.set_xlim((min(x_vals) - margin_pcnt * np.ptp(x_vals), 
                      (max(x_vals) + margin_pcnt * np.ptp(x_vals)))
                    )
        ax.set_ylim((min(y_vals) - margin_pcnt * np.ptp(y_vals), 
                      (max(y_vals) + margin_pcnt * np.ptp(y_vals)))
                    )


    # Animation frame-by-frame function
    def animate_frames(frame):
        # Dynamically modify both x and y data
        x_data = x_vals[:frame + 1]  # Modify x data to include points up to the current frame
        y_data = y_vals[:frame + 1]  # Modify y data accordingly
        line.set_data(x_data, y_data)  # Update the line data
        return line,


    # Create animation
    frames = len(x_vals)
    ani = animation.FuncAnimation(fig, animate_frames, frames=frames, interval=200, blit=True)
    if filename:
        ani.save(filename, dpi=300, writer=PillowWriter(fps= frames / length_secs))
    return ani

#%% Script part
def main():
    #%%% Data Setup
    r = np.arange(start = 40, stop = 400 + 40, step = 40)
    E_states = np.array([3025.0,
                          1583.0,
                          1583.0,
                          3157.0,
                          3157.0,
                          3157.0,
                          1367.0,
                          1367.0,
                          1367.0])
    s = len(E_states)
    conversion = 83.5934788 #kJ/mol to cm^-1 conversion

    # problem 11-11 from kinetics HW7
    SumOfStates = []
    for E in r:
        E = E * conversion
        G = (E)**s / ( math.factorial(s) * np.prod(E_states) )
        SumOfStates.append(G)

    #%%% Plotting
    # set  dark background
    plt.style.use('dark_background')

     # colors :)
    dblue = sns.color_palette("colorblind")[3]

    # Create figure and axis
    fig, ax = plt.subplots()

    # Initial plot setup
    plottyplot = ax.plot(r, SumOfStates, marker='o', linestyle='-', linewidth=4, color=dblue)
    ax.tick_params(axis='x', labelsize=16, pad=8)
    ax.tick_params(axis='y', labelsize=22, pad=8)
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    ax.set_ylabel('$G(E)$', fontsize=24, rotation=0, labelpad=30)
    ax.set_xlabel('$E$', fontsize=24, rotation=0, labelpad=30)
    ax.grid(alpha=0.2, linewidth=2, zorder=0)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout(pad=2)
       

    animate(plottyplot, fig, ax, filename = "animation.gif", margin_pcnt = 0.2, length_secs = 10) #!!!

    # Display animation
    fig.show()


if __name__ == "__main__":
    import math
    import seaborn as sns
    main()
    