from pylab import *
import os
from matplotlib.animation import FuncAnimation, ArtistAnimation

tmpdir = "/tmp/tmp0fdtd"


class Detector:
    def __init__(self, capture_period=8, frame_interval=50):
        self.timing_E = None
        self.timing_H = None
        self.timing_S = None
        self.timing_U = None
        
        self.slicing = None
        
        self.grid = None
        self.index = None
        
        self.capture_period = int(capture_period)
        if self.capture_period <= 0:
            raise ValueError("Capture period has to be positive and greater than 0")
        self.frame_interval = frame_interval
    
    def create_fields_template(self, captures):
        if not os.path.exists(tmpdir):
            os.mkdir(tmpdir)
        TYPE = complex128 if self.grid.is_phasor else float64 #Type for memmaps
        sh = tuple([captures] + list(self.grid.E[self.slicing].shape))
        self.timing_E = memmap(os.path.join(tmpdir, "%i_E.dat" % self.index), dtype=TYPE, mode="w+", shape=sh)
        self.timing_H = memmap(os.path.join(tmpdir, "%i_H.dat" % self.index), dtype=TYPE, mode="w+", shape=sh)
        self.timing_S = memmap(os.path.join(tmpdir, "%i_S.dat" % self.index), dtype=float64, mode="w+", shape=sh[:-1])
        self.timing_U = memmap(os.path.join(tmpdir, "%i_U.dat" % self.index), dtype=float64, mode="w+", shape=sh[:-1])
    
    def save_fields(self, time_step):
        saved_E = copy(self.grid.E)[self.slicing]
        saved_H = copy(self.grid.H)[self.slicing]
        saved_S = sqrt((saved_E[...,1].real*saved_H[...,2].real-saved_E[...,2].real*saved_H[...,1].real)**2
                      +(saved_E[...,2].real*saved_H[...,0].real-saved_E[...,0].real*saved_H[...,2].real)**2
                      +(saved_E[...,0].real*saved_H[...,1].real-saved_E[...,1].real*saved_H[...,0].real)**2)
        saved_U = (abs(saved_E[...,0])**2*self.grid.permittivity[self.slicing][...,0]
                  +abs(saved_E[...,1])**2*self.grid.permittivity[self.slicing][...,1]
                  +abs(saved_E[...,2])**2*self.grid.permittivity[self.slicing][...,2]
                  +abs(saved_H[...,0])**2*self.grid.permeability[self.slicing][...,0]
                  +abs(saved_H[...,1])**2*self.grid.permeability[self.slicing][...,1]
                  +abs(saved_H[...,2])**2*self.grid.permeability[self.slicing][...,2]) * 0.5
        capture = time_step//self.capture_period
        self.timing_E[capture,...] = saved_E
        self.timing_H[capture,...] = saved_H
        self.timing_S[capture,...] = saved_S
        self.timing_U[capture,...] = saved_U
    
    def visualize(self, scalar_field, **kwargs):
        if self.grid.dim == 3 and all([s != 0 or s != slice(None,None,None)]):
            raise AssertionError("Three dimensional detectors cannot be visualized")
            quit()
        dim = len(self.timing_E.shape)-2 #Already squeezed, substract time dimension and dimension added for vectorial character of field
        if dim == 0:
            self.plot_point(scalar_field, **kwargs)
        if dim == 1:
            self.anim_line(scalar_field, **kwargs)
        if dim == 2:
            self.anim_surface(scalar_field, **kwargs)
    
    def visualizeEx(self):
        self.visualize(array(self.timing_E)[...,0].real
                        ,cmap="jet")
    
    def visualizeEy(self):
        self.visualize(array(self.timing_E)[...,1].real
                        ,cmap="jet")
    
    def visualizeEz(self):
        self.visualize(array(self.timing_E)[...,2].real
                        ,cmap="jet")
    
    def visualizeHx(self):
        self.visualize(array(self.timing_H)[...,0].real
                        ,cmap="jet")
    
    def visualizeHy(self):
        self.visualize(array(self.timing_H)[...,1].real
                        ,cmap="jet")
    
    def visualizeHz(self):
        self.visualize(array(self.timing_H)[...,2].real
                        ,cmap="jet")
    
    def visualizePoynting(self):
        self.visualize(array(self.timing_S)
                        ,cmap="hot"
                        ,fmin=0)
    
    def visualizeEnergy(self):
        self.visualize(array(self.timing_U)
                        ,cmap="hot"
                        ,fmin=0)
    
    def plot_point(self, field, **kwargs):
        fig, ax = subplots()
        ax.plot(field, lw=0.75)
        show()
    
    def anim_line(self, field, **kwargs):
        fig, ax = subplots()
        line, = ax.plot(field[0], lw=0.75)
        ax.set_ylim(field.min()*1.05, field.max()*1.05)
        def update(q):
            line.set_ydata(field[q])
            return line,
        ani = FuncAnimation(fig, update, frames=field.shape[0], interval=self.frame_interval)
        show()
    
    def anim_surface(self, field, **kwargs):
        fig, ax = subplots()
        absMax = lambda f: max(abs(f.max()), abs(f.min()))
        art = ax.imshow(field[0], cmap=kwargs["cmap"])#, vmin=-absMax(field[0]), vmax=absMax(field[0]))
        ax.text(0, 1.05, '', transform = ax.transAxes)
        cb = fig.colorbar(art)
        def update(q):
            art.set_data(field[q])
            del ax.texts[0]
            ax.text(0, 1.05, "Time step %i" % (q*self.capture_period), transform = ax.transAxes)
            fmin = -absMax(field[q]) if not "fmin" in kwargs.keys() else kwargs["fmin"]
            fmax = absMax(field[q]) if not "fmax" in kwargs.keys() else kwargs["fmax"]
            art.set_clim(fmin, fmax)
            return art,
        ani = FuncAnimation(fig, update, frames=field.shape[0], interval=self.frame_interval)
        show()
