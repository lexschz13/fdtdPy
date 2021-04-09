from pylab import *


class Object:
    def __init__(self, center, permittivity=1, permeability=1, conductivity=0, monductivity=0, x_axis=None, y_axis=None, z_axis=None):
        self.center = center
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.z_axis = z_axis
        
        self.permittivity = permittivity * array([1,1,1])
        self.permeability = permeability * array([1,1,1])
        self.refractive_index = sqrt(self.permittivity*self.permeability)
        self.conductivity = conductivity * array([1,1,1])
        self.monductivity = monductivity * array([1,1,1])
        
        self.name = None
        
        self.slicing = None
        
        self.grid = None


class Ellipsoid(Object):
    def gen_mask(self):
        mkx,mky,mkz = [zeros(self.grid.x.shape)]*3
        if self.grid.Nx > 1:
            mkx += (self.grid.x - self.center[0])**2 / self.x_axis**2
        if self.grid.Ny > 1:
            mky += (self.grid.y - self.center[1])**2 / self.y_axis**2
        if self.grid.Nx > 1:
            mkz += (self.grid.z - self.center[2])**2 / self.z_axis**2
        self.slicing = where((mkx+mky+mkz)<=1)


class Rectangle(Object):
    def gen_mask(self):
        mkx,mky,mkz = [True]*3
        if self.grid.Nx > 1:
            mkx = (self.grid.x - self.center[0])**2 <= self.x_axis**2
        if self.grid.Ny > 1:
            mky = (self.grid.y - self.center[1])**2 <= self.y_axis**2
        if self.grid.Nx > 1:
            mkz = (self.grid.z - self.center[2])**2 <= self.z_axis**2
        self.slicing = where(mkx*mky*mkz)
