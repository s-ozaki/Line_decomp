from numpy import sqrt, log, exp

light_speed = 2.9979e5 # km/s
fwhm_to_sigma = 0.5 / sqrt(2.0 * log(2.0)) # scale factor from FWHM to sigma

class funcs:
    # If you can add a new function, you must add a scaling rule in the class, 
    # fit_database in dpio.py
    
    def gaussian(self, x, a, b,c):
        c = c / light_speed * b # velocity width to wavelength width
        c = c * fwhm_to_sigma # FWHM to sigma
        yval = a * exp(-0.5*((x-b)/c)**2)
        return(yval)

    def linear(self, x,a,b):
        yval = a * x + b
        return(yval)

    def lorentzian(self, x, a,b,c):
        c = c / light_speed * b # velocity width to wavelength width
        yval = a *c**2/4.0 / ((x-b)**2 + c**2/4.0)
        return(yval)
        
    def scale(self, scale, component_num, funcinfo, params):
        # This function defines the scaling rule for each function listed in this class.
        # Which parameter should be scaled?
        for i in range(component_num):
            if funcinfo[i][0] == 'linear':
                params[i][0][0] = params[i][0][0]/scale
                params[i][0][1] = params[i][0][1]/scale
                params[i][0][2] = params[i][0][2]/scale
                params[i][0][4] = params[i][0][4]/scale

                params[i][1][0] = params[i][1][0]/scale
                params[i][1][1] = params[i][1][1]/scale
                params[i][1][2] = params[i][1][2]/scale
                params[i][1][4] = params[i][1][4]/scale

            elif funcinfo[i][0] == 'gaussian':
                params[i][0][0] = params[i][0][0]/scale
                params[i][0][1] = params[i][0][1]/scale
                params[i][0][2] = params[i][0][2]/scale
                params[i][0][4] = params[i][0][4]/scale

            elif funcinfo[i][0] == 'lorentzian':
                params[i][0][0] = params[i][0][0]/scale
                params[i][0][1] = params[i][0][1]/scale
                params[i][0][2] = params[i][0][2]/scale
                params[i][0][4] = params[i][0][4]/scale
                

        return(params)
        