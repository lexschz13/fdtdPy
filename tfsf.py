from pylab import *
from .vacuum_constants import *



class TFSF:
    def __init__(self, source, distance_to_source):
        self.source = self.check_source(source)
        self.distance = distance_to_source
        self.side = 1 if self.distance>=0 else -1
            
    def check_source(self, source):
        """
        Checks if source has one dimension less than grid and if it is filling all non-flat dimensions
        """
        if source.direction is None:
            raise AssertionError("Source without direction cannot have TFSF boundary")
            quit()
        if source.grid.dim == 1:
            return source
        
        scsh = source.grid.x[source.slicing].shape
        gsh = (source.grid.Nx, source.grid.Ny, source.grid.Nz)
        
        elif source.grid.dim == 2:
            flat_axis_idx = [N==1 for N in gsh].index(True)
            i,j = (flat_axis_idx-1)%3, (flat_axis_idx+1)%3
            N0,N1 = (gsh[i], gsh[j]) if i<j else (gsh[j],gsh[i])
            if (scsh[0]==1 and scsh[1]==N1) or (scsh[0]==N0 and scsh[1]==1):
                return source
            else:
                raise AssertionError("Source is not filling all grid's line")
                quit()
        
        elif source.grid.dim == 3:
            is_same_sh = [scsh[l] == gsh[l] for l in range(3)]
            if is_same_sh.count(False)==1 and scsh[is_same_sh.index(False)]==1:
                return source
            else:
                raise AssertionError("Source is not filling all grid's line")
                quit()
