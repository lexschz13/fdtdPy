from pylab import *
from .vacuum_constants import *



class TFSF:
    def __init__(self, source, distance_to_source):
        self.distance = distance_to_source
        self.side = 1 if self.distance>=0 else -1
        self.source, self.slicing = self.check_source(source)
        self.amplitude = self.source.amplitude / VACUUM_IMPEDANCE
        
        px,py,pz = self.source.polarization
        sx,sy,sz = self.source.direction
        qx = sy*pz - sz*py
        qy = sz*px - sx*pz
        qz = sx*py - sy*px
        self.polarization = [qx,qy,qz] #Polarization for magnetic field
            
    def check_source(self, source):
        """
        Checks if source has one dimension less than grid and if it is filling all non-flat dimensions
        """
        
        if source.direction is None:
            raise AssertionError("Source without direction cannot have TFSF boundary")
            quit()
        
        gsh = (source.grid.Nx, source.grid.Ny, source.grid.Nz)
        
        if source.grid.dim == 1:
            none_flat_axis_idx = [N!=1 for N in gsh].index(True)
            sl = list(source.slicing)
            sl[none_flat_axis_idx] %= gsh[none_flat_axis_idx]
            sl[none_flat_axis_idx] += self.distance
            if sl[none_flat_axis_idx]<0 or sl[none_flat_axis_idx]>=gsh[none_flat_axis_idx]:
                raise AssertionError("TFSF outside of grid")
                quit()
            return source, tuple(sl)
        
        elif source.grid.dim == 2:
            scsh = source.grid.x[source.slicing].shape
            flat_axis_idx = [N==1 for N in gsh].index(True)
            i,j = (flat_axis_idx-1)%3, (flat_axis_idx+1)%3
            N0,N1 = (gsh[i], gsh[j])
            if scsh[0]==N0 or scsh[0]==N1:
                sl = list(source.slicing)
                int_idx = i if scsh[0]==N1 else j
                sl[int_idx] %= gsh[int_idx]
                sl[int_idx] += self.distance
                if sl[int_idx]<0 or sl[int_idx]>=gsh[int_idx]:
                    raise AssertionError("TFSF outside of grid")
                    quit()
                return source, tuple(sl)
            else:
                raise AssertionError("Source is not filling all grid's line")
                quit()
        
        elif source.grid.dim == 3:
            scsh = source.grid.x[source.slicing].shape
            sl = list(source.slicing)
            int_idx = [type(s) is int for s in sl].index(True)
            i,j = (int_idx-1)%3,(int_idx+1)%3
            if (scsh[0]==gsh[i] and scsh[1]==gsh[j]) or (scsh[0]==gsh[j] and scsh[1]==gsh[i]):
                sl[int_idx] %= gsh[int_idx]
                sl[int_idx] += self.direction
                if sl[int_idx]<0 or sl[int_idx] >= gsh[int_idx]:
                    raise AssertionError("TFSF outside of grid")
                    quit()
                return source, tuple(sl)
            else:
                raise AssertionError("Source is not filling all grid's surface")
                quit()
    
    def temp_func(self, time):
        time -= self.side * self.source.grid.refractive_index[self.slicing]/VACUUM_LIGHT_SPEED * (
                        self.source.direction[0]*(self.source.grid.x[self.slicing][...,None] - self.source.grid.x[self.source.slicing][...,None] - 0.5*self.source.grid.spatial_step) +
                        self.source.direction[1]*(self.source.grid.y[self.slicing][...,None] - self.source.grid.y[self.source.slicing][...,None] - 0.5*self.source.grid.spatial_step) +
                        self.source.direction[2]*(self.source.grid.z[self.slicing][...,None] - self.source.grid.z[self.source.slicing][...,None] - 0.5*self.source.grid.spatial_step))
        time -= 0.5*self.source.grid.time_pace
        return self.source.temp_func(time, **self.source.temp_func_params)
    
    def update_H(self, time):
        self.source.grid.H[self.slicing][...,0] -= self.amplitude * self.polarization[0] * self.temp_func(time)[...,0]
        self.source.grid.H[self.slicing][...,1] -= self.amplitude * self.polarization[1] * self.temp_func(time)[...,1]
        self.source.grid.H[self.slicing][...,2] -= self.amplitude * self.polarization[2] * self.temp_func(time)[...,2]
