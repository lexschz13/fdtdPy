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
        
        scsh = source.grid.x[source.slicing].shape
        gsh = (source.grid.Nx, source.grid.Ny, source.grid.Nz)
        
        if source.grid.dim == 1:
            none_flat_axis_idx = [N!=1 for N in gsh].index(True)
            sl = list(source.slicing)
            sl[0] %= gsh[none_flat_axis_idx]
            sl[0] += self.distance
            if sl[0]<0 or sl[0]>=gsh[none_flat_axis_idx]:
                raise AssertionError("TFSF outside of grid")
                quit()
            return source, tuple(sl)
        
        elif source.grid.dim == 2:
            flat_axis_idx = [N==1 for N in gsh].index(True)
            i,j = (flat_axis_idx-1)%3, (flat_axis_idx+1)%3
            N0,N1 = (gsh[i], gsh[j]) if i<j else (gsh[j],gsh[i])
            if (scsh[0]==1 and scsh[1]==N1) or (scsh[0]==N0 and scsh[1]==1):
                sl = list(source.slicing)
                int_idx = [s is int for s in sl].index(True)
                g_int_idx = min(i,j) if int_idx==0 else max(i,j)
                sl[int_idx] %= gsh[g_int_idx]
                sl[int_idx] += self.distance
                if sl[int_idx]<0 or sl[int_idx]>=gsh[g_int_idx]:
                    raise AssertionError("TFSF outside of grid")
                    quit()
                return source, tuple(sl)
            else:
                raise AssertionError("Source is not filling all grid's line")
                quit()
        
        elif source.grid.dim == 3:
            is_same_sh = [scsh[l] == gsh[l] for l in range(3)]
            if is_same_sh.count(False)==1 and scsh[is_same_sh.index(False)]==1:
                sl = list(source.slicing)
                sl[is_same_sh.index(False)] %= gsh[is_same_sh.index(False)]
                sl[is_same_sh.index(False)] += self.direction
                if sl[is_same_sh.index(False)]<0 or sl[is_same_sh.index(False)] >= gsh[is_same_sh.index(False)]:
                    raise AssertionError("TFSF outside of grid")
                    quit()
                return source, tuple(sl)
            else:
                raise AssertionError("Source is not filling all grid's surface")
                quit()
    
    def temp_func(self, time):
        time -= self.side * self.source.grid.refractive_index[self.slicing]/VACUUM_LIGHT_SPEED * (
                        self.direction[0]*(self.grid.x[self.slicing][...,None] - self.grid.x[self.source.slicing][...,None] - 0.5*self.source.grid.spatial_step) +
                        self.direction[1]*(self.grid.y[self.slicing][...,None] - self.grid.y[self.source.slicing][...,None] - 0.5*self.source.grid.spatial_step) +
                        self.direction[2]*(self.grid.z[self.slicing][...,None] - self.grid.z[self.source.slicing][...,None] - 0.5*self.source.grid.spatial_step))
        time -= 0.5*self.source.grid.time_pace
        return self.source.temp_func(time, **self.source.temp_params)
    
    def update_H(self, time):
        self.source.grid.H[self.slicing][...,0] -= self.amplitude * self.polarization[0] * self.temp_func(time)[...,0]
        self.source.grid.H[self.slicing][...,1] -= self.amplitude * self.polarization[0] * self.temp_func(time)[...,1]
        self.source.grid.H[self.slicing][...,2] -= self.amplitude * self.polarization[0] * self.temp_func(time)[...,2]
