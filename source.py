from pylab import *


def gaussian(t, **kwargs):
    s = 2e-7 if not "sigma" in kwargs.keys() else kwargs["sigma"]
    d = 1e-8 if not "delay" in kwargs.keys() else kwargs["delay"]
    x = (t-d)/s
    return exp(-x**2/2)


def ricker(t, **kwargs):
    s = 2e-7 if not "sigma" in kwargs.keys() else kwargs["sigma"]
    d = 1e-8 if not "delay" in kwargs.keys() else kwargs["delay"]
    x = (t-d)/s
    return (1-x**2)*exp(-x**2/2)


def sinusoid(t, **kwargs):
    f = 1e15 if not "frequency" in kwargs.keys() else kwargs["frequency"]
    p = 0 if not "shift" in kwargs.keys() else kwargs["shift"]
    return sin(2*pi*f*t-p)


def singaussian(t, **kwargs):
    s = 2e-7 if not "sigma" in kwargs.keys() else kwargs["sigma"]
    d = 1e-8 if not "delay" in kwargs.keys() else kwargs["delay"]
    f = 1e15 if not "frequency" in kwargs.keys() else kwargs["frequency"]
    p = 0 if not "shift" in kwargs.keys() else kwargs["shift"]
    return gaussian(t, sigma=s, delay=d) * sinusoid(t, frequency=f, shift=p)

def delta(t, **kwargs):
    t0 = 0 if not "center" in kwargs.keys() else kwargs["center"]
    return 1 if t==t0 else 0



class Source:
    def __init__(self, amplitude=1, polarization=[0,0,1], temp_func=gaussian, direction=None, **kwargs):
        self.amplitude = amplitude 
        self.polarization = polarization
        self.temp_func = temp_func
        self.temp_func_params = kwargs
        self.direction = direction
        
        self.slicing = None
        
        self.grid = None
        
    def update_E(self, time):
        self.grid.E[self.slicing][0] += self.polarization[0] * self.amplitude * self.temp_func(time, **self.temp_func_params)
        self.grid.E[self.slicing][1] += self.polarization[1] * self.amplitude * self.temp_func(time, **self.temp_func_params)
        self.grid.E[self.slicing][2] += self.polarization[2] * self.amplitude * self.temp_func(time, **self.temp_func_params)


class GaussianSource(Source):
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, sigma=1e-14, delay=4e-14):
        super().__init__(amplitude, polarization, temp_func=gaussian, direction=direction, sigma=sigma, delay=delay)
    
    @property
    def sigma(self):
        return self.temp_func["sigma"]
    
    @sigma.setter
    def sigma(self, sigma):
        self.temp_func["sigma"] = sigma
    
    @property
    def delay(self):
        return self.temp_func["delay"]
    
    @delay.setter
    def delay(self, delay):
        self.temp_func["delay"] = delay


class RickerSource(Source):
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, sigma=1e-14, delay=4e-14):
        super().__init__(amplitude, polarization, temp_func=ricker, direction=direction, sigma=sigma, delay=delay)
    
    @property
    def sigma(self):
        return self.temp_func["sigma"]
    
    @sigma.setter
    def sigma(self, sigma):
        self.temp_func["sigma"] = sigma
    
    @property
    def delay(self):
        return self.temp_func["delay"]
    
    @delay.setter
    def delay(self, delay):
        self.temp_func["delay"] = delay


class SinusoidSource(Source):
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, frequency=6e14, shift=0):
        super().__init__(amplitude, polarization, temp_func=sinusoid, direction=direction, frequency=frequency, shift=shift)
    
    @property
    def frequency(self):
        return self.temp_func_params["frequency"]
    
    @frequency.setter
    def sigma(self, frequency):
        self.temp_func_params["frequency"] = frequency
    
    @property
    def shift(self):
        return self.temp_func["shift"]
    
    @shift.setter
    def sigma(self, shift):
        self.temp_func["shift"] = shift


class SinGaussianSource(Source):
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, sigma=1e-14, delay=4e-14, frequency=6e14, shift=0):
        super().__init__(amplitude, polarization, temp_func=singaussian, direction=direction, sigma=sigma, delay=delay, frequency=frequency, shift=shift)
    
    @property
    def sigma(self):
        return self.temp_func["sigma"]
    
    @sigma.setter
    def sigma(self, sigma):
        self.temp_func["sigma"] = sigma
    
    @property
    def delay(self):
        return self.temp_func["delay"]
    
    @delay.setter
    def delay(self, delay):
        self.temp_func["delay"] = delay
    
    @property
    def frequency(self):
        return self.temp_func["frequency"]
    
    @frequency.setter
    def sigma(self, frequency):
        self.temp_func["frequency"] = frequency
    
    @property
    def shift(self):
        return self.temp_func["shift"]
    
    @shift.setter
    def sigma(self, shift):
        self.temp_func["shift"] = shift

class DeltaSource(Source):
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, center=0):
        super().__init__(amplitude, polarization, temp_func=delta, center=center)
    
    @property
    def center(self):
        return self.temp_func["center"]
    
    @center.setter
    def center(self, center):
        self.temp_func["center"] = center
