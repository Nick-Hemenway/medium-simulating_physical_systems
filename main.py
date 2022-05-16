import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class ContactEvent():
    """Callable class that returns zero when the ball engages/disengages contact with the ground.
    """
    
    def __init__(self, r, direction=0):
        """
        Parameters
        ----------
        r : float
            Radius of the ball in meters
        direction : int, optional
            Direction of a zero crossing to trigger the contact event.
            Negative for the ball coming into contact with the ground.
            Positive for the ball leaving contact with the ground, by default 0
        """
        
        self.r = r
        self.direction = direction
        self.terminal = True #terminal is True so that simulation will end on contact event
    
    def __call__(self,t,y):
        """Computes the height of the ball above being in contact
        
        Notes
        -----
        The ball will engage/disengage contact when the height of the center
        of the ball equals the radius of the ball.

        Parameters
        ----------
        t : float
            time in the simulation
        y : array-like
            vector of ball's state variables (height and velocity)

        Returns
        -------
        float
            height above being in contact
        """
        #unpack height and velocity of ball
        x1, x2 = y
        return x1 - self.r

class BouncingBall():
    """Class to simulate a bouncing ball

    Methods
    -------
    in_air(t, y) : Computes the ball's state derivatives while in air
    
    in_contact(t, y) : Computes the ball's state derivatives while in contact
    
    simulate(t_span, x0, max_step) : Simulates the ball bouncing
    """
    
    g = 9.81 #gravity constant [m/s^2]
    
    def __init__(self, m, k, c, ball_radius_cm):
        """
        Parameters
        ----------
        m : float
            ball mass in kg
        k : float
            ball stiffness in N/m
        c : float
            ball viscous damping coef. in N/(m/s)
        ball_radius_cm : float
            ball radius in centimeters
        """
        
        self.r = ball_radius_cm/100 #radius in meters
        self.m = m
        self.k = k
        self.c = c
        
        #create event functions for the ball engaging and disengaging contact
        #note that when coming into contact with the ground, the direction of the
        #zero crossing will be negative (height abve the ground is transitioning to negative)
        #whereas when the ball leaves the ground the zero crossing will be positive
        #(height above the ground is transitioning from negative to positive)
        self.hitting_ground = ContactEvent(self.r, direction=-1)
        self.leaving_ground = ContactEvent(self.r, direction=1)
        
    def in_air(self,t,y):
        """computes the ball's state derivatives while in air

        Parameters
        ----------
        t : float
            simulation time
        y : array-like
            vector containing the ball's state variables (height and velocity)

        Returns
        -------
        list
            vector containing the ball's state derivatives
        """
        
        x1, x2 = y
        return [x2, -self.g]
    
    def in_contact(self,t,y):
        """computes the ball's state derivatives while in contact with the ground

        Parameters
        ----------
        t : float
            simulation time
        y : array-like
            vector containing the ball's state variables (height and velocity)

        Returns
        -------
        list
            vector containing the ball's state derivatives
        """
        
        x1, x2 = y
        x1_dot = x2
        x2_dot = -(1/self.m)*(self.c*x2 + self.k*(x1-self.r) + self.m*self.g)
        return [x1_dot, x2_dot]
    
    def simulate(self, t_span, x0, max_step=0.01):
        """Simulates the time evolution of the ball bouncing

        Parameters
        ----------
        t_span : two element tuple/list
            starting and stopping time of the simulation, e.g. (0,10)
        x0 : array-like
            initial conditions of ball (height and velocity in m and m/s)
        max_step : float, optional
            max step size in simulation, by default 0.01

        Returns
        -------
        tuple
            tuple containing the time vector and height of the ball 
        """
        
        #check the initial height of the ball to determine if it's starting in the air
        in_air = x0[0] > self.r 
        #extract out initial starting and stopping times
        t_start, t_stop = t_span
        
        #create lists with initial conditions that we can append the piecewise solutions to
        t_lst = [t_start]
        z_lst = [x0[0]]
        
        #loop until we reach the desired stopping time
        while t_start < t_stop:
            
            # Here we simulate the ball forward in time using either the air or contact model.
            # Each of these subroutines will terminate when either 1) the final desired simulation stop time is reached,
            # or 2) when a contact event is triggered. At a high level, we are simulating our system forward in time using 
            # the relevant physics model (that depends on the state of the system). The simulation will alternate
            # between the "in_air" model and "in_contact" model switching between the two each time a contact event is triggered

            if in_air:
                sol = solve_ivp(self.in_air, [t_start, t_stop], x0, events=[self.hitting_ground], max_step=max_step)
                
            else:
                sol = solve_ivp(self.in_contact, [t_start, t_stop], x0, events=[self.leaving_ground], max_step=max_step) 
                
            #append solution and time array to list of solutions. Note that the starting time of each solution
            #is the stopping time of the previous solution. To avoid having duplicate time points, we will not include the first
            #data point of each simulation. This is also why we created our solution lists above with the initial conditions already
            #in them
            t_lst.append(sol.t[1::])
            z_lst.append(sol.y[0,1::])
            
            #set the starting time and initial conditions to the stopping time and end condtions of the previous loop
            t_start = sol.t[-1]
            x0 = sol.y[:,-1].flatten()
            
            #if we haven't reached the stopping time yet in the current loop, we must be switching between being in the
            #the air and being in contact
            if t_start < t_stop:
                in_air = not in_air

        #concatenate all of the solutions into a single numpy array
        t = np.hstack(t_lst)
        z = np.hstack(z_lst)
        
        return t,z
            
   
b = BouncingBall(1, 10e3, 1, ball_radius_cm=6)

t,z = b.simulate([0,8], [2, 0])

fig, ax = plt.subplots()
ax.plot(t,z)
