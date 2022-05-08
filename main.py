import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class ContactEvent():
    
    def __init__(self, r, direction=0):
        
        self.r = r
        self.direction = direction
        self.terminal = True
    
    def __call__(self, t,y):
        x1, x2 = y
        return x1 - self.r

class BouncingBall():
    
    g = 9.81 #m/s
    
    def __init__(self, m, k, c, ball_radius_cm=6):
        
        self.r = ball_radius_cm/100 #radisu in meters
        self.m = m
        self.k = k
        self.c = c
        
        self.hitting_ground = ContactEvent(self.r, direction=-1)
        self.leaving_ground = ContactEvent(self.r, direction=1)
        
    def in_air(self,t,y):
        
        x1, x2 = y
        return [x2, -self.g]
    
    def in_contact(self,t,y):
        
        x1, x2 = y
        x1_dot = x2
        x2_dot = -1/self.m*(self.c*x2 + self.k*x1)
        return [x1_dot, x2_dot]
    
    def simulate(self, t_span, x0, max_step=0.01):
        
        in_air = True
        t_start, t_stop = t_span
        t_lst = []
        z_lst = []
        
        while t_start < t_stop:
            
            #simulate ball forward in time using either air or contact model
            if in_air:
                sol = solve_ivp(self.in_air, [t_start, t_stop], x0, events=[self.hitting_ground], max_step=max_step)
                
            else:
                sol = solve_ivp(self.in_contact, [t_start, t_stop], x0, events=[self.leaving_ground], max_step=max_step) 
                
            #append solution and time array to list of solutions
            t_lst.append(sol.t)
            z_lst.append(sol.y[0,:])
            
            #reset the initial conditions and starting time
            t_start = sol.t[-1]
            x0 = sol.y[:,-1].flatten()
            
            if t_start < t_stop:
                in_air = not in_air

        t = np.hstack(t_lst)
        z = np.hstack(z_lst)
        
        return t,z
            
   
b = BouncingBall(1, 10e3, 1)

t,z = b.simulate([0,8], [2, 0])

fig, ax = plt.subplots()
ax.plot(t,z)
    
    
    