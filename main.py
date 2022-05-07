import numpy as np
from scipy.integrate import solve_ivp

class BouncingBall():
    
    g = 9.81 #m/s
    
    def __init__(self, m, k, c, ball_radius_cm=6):
        
        self.r = ball_radius_cm/100 #radisu in meters
        self.m = m
        self.k = k
        self.c = c
        
    def in_air(self,t,y):
        
        x1, x2 = y
        return [x2, -BouncingBall.g]
    
    def in_contact(self,t,y):
        
        x1, x2 = y
        x1_dot = x2
        x2_dot = -1/self.m*(self.c*x2 + self.k*x1)
        return [x1_dot, x2_dot]
    
    
    