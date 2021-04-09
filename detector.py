from pylab import *
from matplotlib.animation import FuncAnimation, ArtistAnimation



class Detector:
    def __init__(self):
        self.timing_E = []
        self.timing_H = []
        self.timing_S = []
        self.timing_U = []
        
        self.slicing = None
        
        self.grid = None
    
    def save_fields(self):
        saved_E = copy(self.grid.E)[self.slicing]
        saved_H = copy(self.grid.H)[self.slicing]
        saved_S = sqrt((saved_E[...,1]*saved_H[...,2]-saved_E[...,2]*saved_H[...,1])**2
                       +(saved_E[...,2]*saved_H[...,0]-saved_E[...,0]*saved_H[...,2])**2
                       +(saved_E[...,0]*saved_H[...,1]-saved_E[...,1]*saved_H[...,0])**2)
        saved_U = (saved_E[...,0]**2*self.grid.permittivity[self.slicing][...,0]
                  +saved_E[...,1]**2*self.grid.permittivity[self.slicing][...,1]
                  +saved_E[...,2]**2*self.grid.permittivity[self.slicing][...,2]
                  +saved_H[...,0]**2*self.grid.permeability[self.slicing][...,0]
                  +saved_H[...,1]**2*self.grid.permeability[self.slicing][...,1]
                  +saved_H[...,2]**2*self.grid.permeability[self.slicing][...,2]) * 0.5
        self.timing_E.append(saved_E)
        self.timing_H.append(saved_H)
        self.timing_S.append(saved_S)
        self.timing_U.append(saved_U)
    
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
        self.visualize(array(self.timing_E)[...,0]
                        ,cmap="jet")
    
    def visualizeEy(self):
        self.visualize(array(self.timing_E)[...,1]
                        ,cmap="jet")
    
    def visualizeEz(self):
        self.visualize(array(self.timing_E)[...,2]
                        ,cmap="jet")
    
    def visualizeHx(self):
        self.visualize(array(self.timing_H)[...,0]
                        ,cmap="jet")
    
    def visualizeHy(self):
        self.visualize(array(self.timing_H)[...,1]
                        ,cmap="jet")
    
    def visualizeHz(self):
        self.visualize(array(self.timing_H)[...,2]
                        ,cmap="jet")
    
    def visualizePoynting(self):
        """
        Ex = array(self.timing_E)[...,0]
        Ey = array(self.timing_E)[...,1]
        Ez = array(self.timing_E)[...,2]
        Hx = array(self.timing_H)[...,0]
        Hy = array(self.timing_H)[...,1]
        Hz = array(self.timing_H)[...,2]
        Sx = Ey*Hz - Ez*Hy
        Sy = Ez*Hx - Ex*Hz
        Sz = Ex*Hy - Ey*Hx
        """
        self.visualize(array(self.timing_S)
                        ,cmap="hot"
                        ,fmin=0)
    
    def visualizeEnergy(self):
        """
        Ex = array(self.timing_H)[...,0]
        Ey = array(self.timing_H)[...,1]
        Ez = array(self.timing_H)[...,2]
        Hx = array(self.timing_H)[...,0]
        Hy = array(self.timing_H)[...,1]
        Hz = array(self.timing_H)[...,2]
        """
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
        ani = FuncAnimation(fig, update, frames=field.shape[0], interval=50)
        show()
    
    def anim_surface(self, field, **kwargs):
        fig, ax = subplots()
        absMax = lambda f: max(abs(f.max()), abs(f.min()))
        art = ax.imshow(field[0], cmap=kwargs["cmap"])#, vmin=-absMax(field[0]), vmax=absMax(field[0]))
        cb = fig.colorbar(art)
        def update(q):
            art.set_data(field[q])
            fmin = -absMax(field[q]) if not "fmin" in kwargs.keys() else kwargs["fmin"]
            fmax = absMax(field[q]) if not "fmax" in kwargs.keys() else kwargs["fmax"]
            art.set_clim(fmin, fmax)
            return art,
        ani = FuncAnimation(fig, update, frames=field.shape[0], interval=50)
        show()
