from pylab import *
from .vacuum_constants import *


possible_disp = ["drude"]


class Object:
    def __init__(self, center, permittivity=1, permeability=1, conductivity=0, monductivity=0, x_axis=None, y_axis=None, z_axis=None, dispersion=None, **disp_kwargs):
        self.center = center
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.z_axis = z_axis
        
        self.permittivity = permittivity * array([1,1,1])
        self.permeability = permeability * array([1,1,1])
        self.refractive_index = (abs(self.permittivity*self.permeability))**0.5
        self.conductivity = conductivity * array([1,1,1])
        self.monductivity = monductivity * array([1,1,1])
        
        self.name = None
        
        self.slicing = None
        
        self.grid = None
        
        self.dispersion = dispersion
        self.disp_kwargs = disp_kwargs
    
    def dispersion_coefs(self):
        if self.dispersion is None:
            pass
        
        elif self.dispersion == "drude":
            g = 5 if not "damping" in self.disp_kwargs.keys() else self.disp_kwargs["damping"]
            fp = 1e15 if not "plasma_frequency" in self.disp_kwargs.keys() else self.disp_kwargs["plasma_frequency"]
            eps_inf = self.permittivity
            self.prevE = self.grid.E[self.slicing]
            self.dJp = zeros(self.prevE.shape, dtype=self.prevE.dtype)
            Ng = 1/(g*self.grid.time_pace)
            Np = VACUUM_LIGHT_SPEED/(self.grid.spatial_step*fp)
            self.cjj = (1 + 0.5/Ng)/(1 - 0.5/Ng)
            self.cje = 1/(1-0.5/Ng) * 2*pi**2*self.grid.courant_number / (VACUUM_IMPEDANCE*Np**2)
            
            cee = (eps_inf - 0.5*self.conductivity*self.grid.time_pace/VACUUM_PERMITTIVITY - 0.5*self.cje*VACUUM_IMPEDANCE*self.grid.courant_number) / (eps_inf + 0.5*self.conductivity*self.grid.time_pace/VACUUM_PERMITTIVITY + 0.5*self.cje*VACUUM_IMPEDANCE*self.grid.courant_number)
            self.grid.cexe[self.slicing] = cee[...,0]
            self.grid.ceye[self.slicing] = cee[...,1]
            self.grid.ceze[self.slicing] = cee[...,2]
            
            ceh = self.grid.courant_number*VACUUM_IMPEDANCE / (eps_inf + 0.5*self.conductivity*self.grid.time_pace/VACUUM_PERMITTIVITY + 0.5*self.cje*VACUUM_IMPEDANCE*self.grid.courant_number)
            self.grid.cexh[self.slicing] = ceh[...,0]
            self.grid.ceyh[self.slicing] = ceh[...,1]
            self.grid.cezh[self.slicing] = ceh[...,2]
        
        else:
            raise ValueError("Only the following types of dispersion are allowed: %s" % ", ".join(possible_disp))
            quit()
    
    def update_J(self):
        if self.dispersion is None:
            pass
        
        elif self.dispersion == "drude":
            self.dJp *= self.cjj
            self.dJp += self.cje * (self.grid.E[self.slicing] + self.prevE)
        
        else:
            raise ValueError("Only the following types of dispersion are allowed: %s" % ", ".join(possible_disp))
            quit()


class Ellipsoid(Object):
    def gen_mask(self):
        mkx,mky,mkz = [zeros(self.grid.x.shape)]*3
        if self.grid.Nx > 1 and self.x_axis != None:
            mkx += (self.grid.x - self.center[0])**2 / self.x_axis**2
        if self.grid.Ny > 1 and self.y_axis != None:
            mky += (self.grid.y - self.center[1])**2 / self.y_axis**2
        if self.grid.Nz > 1 and self.z_axis != None:
            mkz += (self.grid.z - self.center[2])**2 / self.z_axis**2
        self.slicing = where((mkx+mky+mkz)<=1)


class Quadrilateral(Object):
    def gen_mask(self):
        mkx,mky,mkz = [True]*3
        if self.grid.Nx > 1 and self.x_axis != None:
            mkx = (self.grid.x - self.center[0])**2 <= self.x_axis**2
        if self.grid.Ny > 1 and self.y_axis != None:
            mky = (self.grid.y - self.center[1])**2 <= self.y_axis**2
        if self.grid.Nz > 1 and self.z_axis != None:
            mkz = (self.grid.z - self.center[2])**2 <= self.z_axis**2
        self.slicing = where(mkx*mky*mkz)
