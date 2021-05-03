from pylab import *
import os
from tempfile import mkdtemp
from tqdm import tqdm

from .vacuum_constants import *
from .detector import Detector, tmpdir
from .boundaries import Boundary
from .sources import Source
from .objects import Object


"""
nbits = 32

if nbits==16:
    C_DATA = copmlex32
    F_DATA = float16
elif nbits==32:
    C_DATA = complex64
    F_DATA = float32
elif nbits = 64:
    C_DATA = complex128
    F_DATA = float64
"""



class Grid:
    def __init__(self, shape, spatial_step = 1.55e-7, courant_number=None, permittivity=1, permeability=1, conductivity=0, monductivity=0, is_phasor=True, coord_centered=True):
        self.Nx,self.Ny,self.Nz = shape
        if self.Nx<1 or self.Ny<1 or self.Nz<1:
            raise ValueError("Shape must be 1 or more")
            quit()
        if not type(self.Nx) is int and type(self.Ny) is int and type(self.Nz) is int:
            raise TypeError("Shape must be expressed in integers")
            quit()
        self.dim = int(self.Nx>1) + int(self.Ny>1) + int(self.Nz>1)
        if self.dim==0:
            raise ValueError("Grid has to be at least one dimension")
            quit()
        self.spatial_step = spatial_step
        self.courant_number = courant_number if courant_number!=None else 1/sqrt(self.dim)
        self.time_pace = self.courant_number * self.spatial_step / VACUUM_LIGHT_SPEED
        
        a,b,c = spatial_step * array([self.Nx,self.Ny,self.Nz])
        if coord_centered:
            self.x,self.y,self.z = meshgrid(linspace(-a/2,a/2-spatial_step,self.Nx),
                                            linspace(-b/2,b/2-spatial_step,self.Ny),
                                            linspace(-c/2,c/2-spatial_step,self.Nz))
        else:
            self.x,self.y,self.z = meshgrid(linspace(0,a-spatial_step,self.Nx),
                                            linspace(0,b-spatial_step,self.Ny),
                                            linspace(0,c-spatial_step,self.Nz))
        self.x = self.x.transpose(1,0,2)
        self.y = self.y.transpose(1,0,2)
        self.z = self.z.transpose(1,0,2)
        
        self.permittivity = ones((self.Nx,self.Ny,self.Nz,3)) * (permittivity * array([1,1,1]))
        self.permeability = ones((self.Nx,self.Ny,self.Nz,3)) * (permeability * array([1,1,1]))
        self.refractive_index = sqrt(abs(self.permittivity*self.permeability))
        self.conductivity = ones((self.Nx,self.Ny,self.Nz,3)) * (conductivity * array([1,1,1]))
        self.monductivity = ones((self.Nx,self.Ny,self.Nz,3)) * (monductivity * array([1,1,1]))
        
        self.chxe, self.chye, self.chze, self.cexh, self.ceyh, self.cezh = [None]*6
        
        self.is_phasor = is_phasor
        F_TYPE = complex128 if self.is_phasor else float64
        self.E = zeros((self.Nx,self.Ny,self.Nz,3), dtype=F_TYPE)
        self.H = zeros((self.Nx,self.Ny,self.Nz,3), dtype=F_TYPE)
        
        self.objects = {}
        self.sources = []
        self.boundaries = {"xH":None, "xE":None, "yH":None, "yE":None, "zH":None, "zE":None}
        self.detectors = []
        
    def __setitem__(self, key, attr):
        if type(key) != str:
            none_sl = slice(None,None,None)
            if not type(key) is tuple:
                key = (key,)
            if len(key) > 3:
                raise AssertionError("Grid lives in a three-dimesional space")
                quit()
            if not all([type(s) is int or type(s) is slice for s in key]):
                raise TypeError("Slice must be an integer or a sequence of integers made by PySlice")
                quit()
            key = self.fill_slicing(key)
            non_flat_sl_assumed = 3-(key.count(none_sl)+key.count(0))
            if self.dim == 1 and non_flat_sl_assumed > 1:
                raise AssertionError("Invalid slicing for one dimension grid")
                quit()
            if self.dim == 2 and non_flat_sl_assumed > 2:
                raise AssertionError("Invalid slicing for two dimension grid")
                quit()
            if self.dim == 3:
                pass
            self.check_slicing(key)
        if issubclass(type(attr), Source):
            self.add_source(attr, key)
        elif issubclass(type(attr), Detector):
            self.add_detector(attr, key)
        elif issubclass(type(attr), Object) and type(key) is str:
            self.add_object(attr, key)
        elif issubclass(type(attr), Boundary):
            self.add_boundary(attr, key)
        else:
            raise TypeError("Wrong type insertion on grid")
    
    def fill_slicing(self, slicing):
        is_flat_axis = [n == 1 for n in [self.Nx, self.Ny, self.Nz]]
        none_sl = slice(None,None,None)
        l = len(slicing)
        new_sl = [none_sl]*3
        if l == 1:
            new_sl[is_flat_axis.index(False)] = slicing[0]
        if l == 2:
            new_sl[(is_flat_axis.index(True)-1)%3] = slicing[0]
            new_sl[(is_flat_axis.index(True)+1)%3] = slicing[1]
        if l == 3:
            new_sl = [s for s in slicing]
        for i in range(3):
            if new_sl[i]==none_sl and is_flat_axis[i]:
                new_sl[i] = 0
        return tuple(new_sl)
    
    def check_slicing(self, slicing):
        Ns = [n for n in [self.Nx,self.Ny,self.Nz]]
        for i,sl in enumerate(slicing):
            if type(sl) is int:
                sl += Ns[i]*(sl<0)
                sl = slice(sl,sl+1,None)
            if not type(sl) is slice:
                raise TypeError("Only integers for slicing")
            start,end,_ = sl.indices(Ns[i])
            start = start + (Ns[i]*(start<0)) if start!=None else 0
            end = end + (Ns[i]*(end<0)) if end!=None else Ns[i]
            if start>=end:
                raise ValueError("Start of slice must be lower than end")
                quit()
            if end>Ns[i] or start<0:
                raise IndexError("Indices out of range")
                quit()
        
    def add_object(self, obj, name):
        obj.grid = self
        obj.gen_mask()
        self.objects[name] = obj
    
    def add_boundary(self, bound, slicing):
        bound.grid = self
        bound.check_slicing(slicing)
        bound.other_coefs()
        self.boundaries[bound.orientation+bound.field] = bound
    
    def add_source(self, source, slicing):
        source.grid = self
        source.slicing = slicing
        source.add_params()
        source.check_tfsf()
        self.sources += [source]
    
    def add_detector(self, detector, slicing):
        detector.grid = self
        detector.slicing = slicing
        detector.index = len(self.detectors)
        self.detectors += [detector]
    
    def set_coefs(self):
        chc = (0.5*self.monductivity/self.permeability)*self.courant_number*self.spatial_step/VACUUM_IMPEDANCE
        che = self.courant_number/self.permeability/VACUUM_IMPEDANCE
        cec = (0.5*self.conductivity/self.permittivity)*self.courant_number*self.spatial_step*VACUUM_IMPEDANCE
        ceh = self.courant_number/self.permittivity*VACUUM_IMPEDANCE
        
        self.chxe = che[...,0]/(1+chc[...,0])
        self.chxh = (1-chc[...,0])/(1+chc[...,0])
        self.chye = che[...,1]/(1+chc[...,1])
        self.chyh = (1-chc[...,1])/(1+chc[...,1])
        self.chze = che[...,2]/(1+chc[...,2])
        self.chzh = (1-chc[...,2])/(1+chc[...,2])
        
        self.cexe = (1-cec[...,0])/(1+cec[...,0])
        self.cexh = ceh[...,0]/(1+cec[...,0])
        self.ceye = (1-cec[...,1])/(1+cec[...,1])
        self.ceyh = ceh[...,1]/(1+cec[...,1])
        self.ceze = (1-cec[...,2])/(1+cec[...,2])
        self.cezh = ceh[...,2]/(1+cec[...,2])
        
        for name in self.objects.keys():
            obj = self.objects[name]
            obj.dispersion_coefs()
    
    def update_H(self):
        #Lossy factor
        self.H[...,0] *= self.chxh
        self.H[...,1] *= self.chyh
        self.H[...,2] *= self.chzh
        
        #Update on x axis
        if self.Nx > 1:
            self.H[:-1,:,:,1] += self.chye[:-1,:,:] * (self.E[1:,:,:,2] - self.E[:-1,:,:,2]) #Hy
            self.H[:-1,:,:,2] -= self.chze[:-1,:,:] * (self.E[1:,:,:,1] - self.E[:-1,:,:,1]) #Hz
        
        #Update on y axis
        if self.Ny > 1:
            self.H[:,:-1,:,0] -= self.chxe[:,:-1,:] * (self.E[:,1:,:,2] - self.E[:,:-1,:,2]) #Hx
            self.H[:,:-1,:,2] += self.chze[:,:-1,:] * (self.E[:,1:,:,0] - self.E[:,:-1,:,0]) #Hz
        
        #Update on z axis
        if self.Nz > 1:
            self.H[:,:,:-1,0] += self.chxe[:,:,:-1] * (self.E[:,:,1:,1] - self.E[:,:,:-1,1]) #Hx
            self.H[:,:,:-1,1] -= self.chye[:,:,:-1] * (self.E[:,:,1:,0] - self.E[:,:,:-1,0]) #Hy
        
        #Update H in boundaries subgrid (PML)
        for bound_key in self.boundaries.keys():
            if not self.boundaries[bound_key] is None:
                self.boundaries[bound_key].update_boundH()
        
        #Apply H boundaries
        if not self.boundaries["xH"] is None:
            self.boundaries["xH"].update_H()
        if not self.boundaries["yH"] is None:
            self.boundaries["yH"].update_H()
        if not self.boundaries["zH"] is None:
            self.boundaries["zH"].update_H()
        
        #TFSF
        if len(self.sources) == 1:
            src = self.sources[0]
            if src.tfsf_for != None:
                src.tfsf_for.update_H(q * self.time_pace)
            if src.tfsf_back != None:
                src.tfsf_back.update_H(q * self.time_pace)
    
    def update_E(self, time_step):
        #Lossy factor
        self.E[...,0] *= self.cexe
        self.E[...,1] *= self.ceye
        self.E[...,2] *= self.ceze
        
        #Update on x axis
        if self.Nx > 1:
            self.E[1:,:,:,1] -= self.ceyh[1:,:,:] * (self.H[1:,:,:,2] - self.H[:-1,:,:,2]) #Ey
            self.E[1:,:,:,2] += self.cezh[1:,:,:] * (self.H[1:,:,:,1] - self.H[:-1,:,:,1]) #Ez
        
        #Update on y axis
        if self.Ny > 1:
            self.E[:,1:,:,0] += self.cexh[:,1:,:] * (self.H[:,1:,:,2] - self.H[:,:-1,:,2]) #Ex
            self.E[:,1:,:,2] -= self.cezh[:,1:,:] * (self.H[:,1:,:,0] - self.H[:,:-1,:,0]) #Ez
        
        #Update on z axis
        if self.Nz > 1:
            self.E[:,:,1:,0] -= self.cexh[:,:,1:] * (self.H[:,:,1:,1] - self.H[:,:,:-1,1]) #Ex
            self.E[:,:,1:,1] += self.ceyh[:,:,1:] * (self.H[:,:,1:,0] - self.H[:,:,:-1,0]) #Ey
        
        for bound_key in self.boundaries.keys():
            if not self.boundaries[bound_key] is None:
                self.boundaries[bound_key].update_boundE()
        
        #Apply E boundaries
        if not self.boundaries["xE"] is None:
            self.boundaries["xE"].update_E()
        if not self.boundaries["yE"] is None:
            self.boundaries["yE"].update_E()
        if not self.boundaries["zE"] is None:
            self.boundaries["zE"].update_E()
        
        #Injection
        for source in self.sources:
            source.update_E(time_step * self.time_pace)
        
        #Disperssion correction
        for name in self.objects:
            obj = self.objects[name]
            if not obj.dispersion is None:
                self.E[obj.slicing][...,0] -= self.cexh[obj.slicing] * 0.5*(1+obj.cjj) * obj.dJp[...,0]
                self.E[obj.slicing][...,1] -= self.ceyh[obj.slicing] * 0.5*(1+obj.cjj) * obj.dJp[...,1]
                self.E[obj.slicing][...,2] -= self.cezh[obj.slicing] * 0.5*(1+obj.cjj) * obj.dJp[...,2]
        
        #Updating Jp in dispersive objects
        for name in self.objects.keys():
            obj = self.objects[name]
            obj.update_J()
        
    def run(self, time_steps):
        if os.path.exists(tmpdir):
            for fl in os.listdir(tmpdir): #Remove past simulations
                os.remove(os.path.join(tmpdir, fl))
        for d in self.detectors:
            d.create_fields_template(int(time_steps/d.capture_period))
        
        for name in self.objects:
            obj = self.objects[name]
            self.permittivity[obj.slicing] = obj.permittivity
            self.permeability[obj.slicing] = obj.permeability
            self.refractive_index[obj.slicing] = obj.refractive_index
            self.conductivity[obj.slicing] = obj.conductivity
            self.monductivity[obj.slicing] = obj.monductivity
        
        self.set_coefs()
            
        for q in tqdm(range(time_steps)):
            self.update_H()
            
            #Saving pervious E for dispersive objects
            for name in self.objects.keys():
                obj = self.objects[name]
                if not obj.dispersion is None:
                    obj.prevE = self.E[obj.slicing]
            
            self.update_E(q)
            
            for d in self.detectors:
                if q%d.capture_period == 0:
                    d.save_fields(q)
                
    def close(self):
        for d in self.detectors:
            d.timing_E._mmap.close()
            d.timing_H._mmap.close()
            d.timing_S._mmap.close()
            d.timing_U._mmap.close()
    
    def visualizeEx(self):
        for d in self.detectors:
            d.visualizeEx()
    
    def visualizeEy(self):
        for d in self.detectors:
            d.visualizeEy()
    
    def visualizeEz(self):
        for d in self.detectors:
            d.visualizeEz()
    
    def visualizeHx(self):
        for d in self.detectors:
            d.visualizeHx()
    
    def visualizeHy(self):
        for d in self.detectors:
            d.visualizeHy()
    
    def visualizeHz(self):
        for d in self.detectors:
            d.visualizeHz()
    
    def visualizePoynting(self):
        for d in self.detectors:
            d.visualizePoynting()
    
    def visualizeEnergy(self):
        for d in self.detectors:
            d.visualizeEnergy()
