import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation

plt.style.use('bmh')

df = pd.read_csv('sim_data.csv')

# first convert data to uniformly sampled at standard fps
fps = 30
t_new = np.arange(0, df['time'].max(), 1/fps)
z_new = np.interp(t_new, df['time'], df['z'])

df = pd.DataFrame({'time':t_new, 'z':z_new}).set_index('time')

class AnimatedBall():
    
    def __init__(self, t, z):
        self.t = t
        self.z = z
        self.fig, self.ax = plt.subplots()
        self.xlim = (-0.02, self.t.max()*1.02)
        self.ylim = (-0.02, self.z.max()*1.02)
        
    def plot_frame(self, frame_number):
        
        self.ax.cla()
        self.ax.set_xlim(*self.xlim)
        self.ax.set_ylim(*self.ylim)
        self.ax.plot(self.t[0:frame_number], self.z[0:frame_number])
        self.ax.plot(self.t[frame_number-1], self.z[frame_number-1], marker='o', ms=10)
        # self.ax.grid()
        self.ax.set_xlabel('Time [s]')
        self.ax.set_ylabel('Height [m]')
        self.fig.tight_layout()
        
        return self.fig, self.ax

an = AnimatedBall(t_new, z_new)

an.plot_frame(150)
anim = FuncAnimation(an.fig, an.plot_frame, frames=len(an.t))
anim.save('bouncing_ball.gif', fps=fps)        
        