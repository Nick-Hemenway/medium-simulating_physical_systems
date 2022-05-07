import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class BouncingBall():
    
    g = 9.81 #m/s
    
    def __init__(self, m, k, c, ball_radius_cm=6):
        
        self.r = ball_radius_cm/100 #radisu in meters
        self.m = m
        self.k = k
        self.c = c
        
        BouncingBall.contact_event.terminal = True
        
    def in_air(self,t,y):
        
        x1, x2 = y
        return [x2, -BouncingBall.g]
    
    def in_contact(self,t,y):
        
        x1, x2 = y
        x1_dot = x2
        x2_dot = -1/self.m*(self.c*x2 + self.k*x1)
        return [x1_dot, x2_dot]
    
    def contact_event(self, t, y):
        x1, x2 = y
        return x1 - self.r
    
    def simulate(self, t_span, x0):
        
        in_air = True
        t_start, t_stop = t_span
        t_lst = []
        z_lst = []
        
        while t_start <= t_stop:
            
            if in_air:
                sol = solve_ivp(self.in_air, [t_start, t_stop], x0, events=[self.contact_event], max_step=0.1)
                in_air = False
            else:
                sol = solve_ivp(self.in_contact, [t_start, t_stop], x0, events=[self.contact_event], max_step=0.1) 
                in_air = True
                
            t_lst.append(sol.t)
            z_lst.append(sol.y[0,:])
            t_start = sol.t[-1]
   
b = BouncingBall(1, 1, 1)

sol = b.simulate([0,10], [2, 0])

fig, ax = plt.subplots()
ax.plot(sol.t, sol.y[0,:])
    
    
    