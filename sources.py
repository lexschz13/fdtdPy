from pylab import *
from .vacuum_constants import *
from .tfsf import TFSF


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


def saturation(t, **kwargs):
    w = 1 if not "width" in kwargs.keys() else kwargs["width"]
    p = 1 if not "pace" in kwargs.keys() else kwargs["pace"]
    if hasattr(t, "__iter__"):
        result = zeros(t.shape)
        wh_not_null = where(t!=0)
        result[wh_not_null] = (p+1)/(p + exp(w/t[wh_not_null]))
    else:
        result = 0 if t==0 else (p+1)/(p + exp(w/t))
    return result


def sinusoid(t, **kwargs):
    f = 1e15 if not "frequency" in kwargs.keys() else kwargs["frequency"]
    p = 0 if not "shift" in kwargs.keys() else kwargs["shift"]
    c = False if not "is_phasor" in kwargs.keys() else kwargs["is_phasor"]
    if c:
        return -1j*exp(1j*(2*pi*f*t-p))
    else:
        return sin(2*pi*f*t-p)



class Source:
    def __init__(self, amplitude=1, polarization=[0,0,1], env_func=gaussian, direction=None, frequency=0, shift=-pi/2, tfsf_for_d=None, tfsf_back_d=None, **kwargs):
        self.amplitude = amplitude 
        self.polarization = polarization
        self.envelope_func = env_func
        self.temp_func_params = kwargs
        self.temp_func_params.update({"frequency":frequency, "shift":shift})
        self.direction = direction
        
        if not tfsf_for_d is None:
            if tfsf_for_d < 0:
                raise ValueError("Distance for forward TFSF has to be positive or 0")
                quite()
        self.tfsf_for_d  = tfsf_for_d
        self.tfsf_for = None
        if not tfsf_back_d is None:
            if tfsf_back_d >= 0:
                raise ValueError("Distance for forward TFSF has to be negative")
                quite()
        self.tfsf_back_d = tfsf_back_d
        self.tfsf_back = None
        
        self.slicing = None
        
        self.grid = None
        
    def update_E(self, time):
        self.grid.E[self.slicing][...,0] += self.polarization[0] * self.amplitude * self.temp_func(time, **self.temp_func_params)[...,0]
        self.grid.E[self.slicing][...,1] += self.polarization[1] * self.amplitude * self.temp_func(time, **self.temp_func_params)[...,1]
        self.grid.E[self.slicing][...,2] += self.polarization[2] * self.amplitude * self.temp_func(time, **self.temp_func_params)[...,2]
    
    def correct_time(self, time):
        if not self.direction is None:
            return time - self.grid.refractive_index[self.slicing]/VACUUM_LIGHT_SPEED * (
                            self.direction[0]*self.grid.x[self.slicing][...,None] +
                            self.direction[1]*self.grid.y[self.slicing][...,None] +
                            self.direction[2]*self.grid.z[self.slicing][...,None])
        else:
            return time * array([1,1,1])
    
    def check_tfsf(self):
        if len(self.grid.sources) >= 1:
            if self.tfsf_for_d != None or self.tfsf_back_d != None:
                raise AssertionError("Only one source if TFSF is used")
                quit()
            if self.grid.sources[0].tfsf_for!=None or self.grid.sources[0].tfsf_back!=None:
                raise AssertionError("Only one source if TFSF is used")
                quit()
        self.tfsf_for = TFSF(self, self.tfsf_for_d) if not self.tfsf_for_d is None else None
        self.tfsf_back = TFSF(self, self.tfsf_back_d) if not self.tfsf_back_d is None else None
    
    def add_params(self):
        self.temp_func_params.update({"is_phasor":self.grid.is_phasor})
    
    def temp_func(self, time, **kwargs):
        return self.envelope_func(time, **kwargs) * sinusoid(self.correct_time(time), **kwargs)
    
    def model_source(self, time, **kwargs):
        return self.envelope_func(time, **kwargs) * sinusoid(time, **kwargs)
    
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


class GaussianSource(Source):
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, frequency=0, shift=-pi/2, sigma=1e-14, delay=4e-14, tfsf_for_d=None, tfsf_back_d=None):
        super().__init__(amplitude, polarization, env_func=gaussian, direction=direction, frequency=frequency, shift=shift, tfsf_for_d=tfsf_for_d, tfsf_back_d=tfsf_back_d, sigma=sigma, delay=delay)
    
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
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, frequency=0, shift=-pi/2, sigma=1e-14, delay=4e-14, tfsf_for_d=None, tfsf_back_d=None):
        super().__init__(amplitude, polarization, env_func=ricker, direction=direction, frequency=frequency, shift=shift, tfsf_for_d=tfsf_for_d, tfsf_back_d=tfsf_back_d, sigma=sigma, delay=delay)
    
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


class SatSource(Source):
    def __init__(self, amplitude=1, polarization=[0,0,1], direction=None, frequency=0, shift=-pi/2, width=1, pace=1, tfsf_for_d=None, tfsf_back_d=None):
        super().__init__(amplitude, polarization, env_func=saturation, frequency=frequency, shift=shift, tfsf_for_d=tfsf_for_d, tfsf_back_d=tfsf_back_d, direction=direction, width=width, pace=pace)
    
    @property
    def width(self):
        return self.temp_func_params["width"]
    
    @width.setter
    def sigma(self, width):
        self.temp_func_params["width"] = width
    
    @property
    def pace(self):
        return self.temp_func["pace"]
    
    @pace.setter
    def sigma(self, pace):
        self.temp_func["pace"] = pace
