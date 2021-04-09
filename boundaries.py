from pylab import *


class Boundary:
    def __init__(self):
        self.grid = None
        
        self.slicing = None
        self.orientation = ""
        self.field = ""
    
    def check_slicing(self, slicing): #Checks slicing corresponds to boundary
        none_sl = slice(None,None,None)
        is_flat_axis = [n ==1 for n in [self.grid.Nx, self.grid.Ny, self.grid.Nz]]
        if self.grid.dim == 1:
            i = is_flat_axis.index(False)
            if not type(slicing[i]) is int:
                raise AssertionError("In one dimension boundary is defined by one point")
                quit()
            if slicing[i] != 0 and slicing[i] != -1:
                raise ValueError("Boundary must be placed at limit of the grid, indices 0 or -1")
                quit()
            else:
                self.slicing = slicing
                self.field = "E" if slicing[i] == 0 else "H"
                self.orient()
        elif self.grid.dim == 2:
            j = (is_flat_axis.index(True)-1)%3
            k = (is_flat_axis.index(True)+1)%3
            if not (type(slicing[j]) is int or type(slicing[k]) is int) and type(slicing[j]) != type(slicing[k]):
                raise AssertionError("In two dimensions boundary is defined by a line")
                quit()
            point_line = slicing[j] if type(slicing[j]) is int else slicing[k]
            if point_line != 0 and point_line != -1:
                raise ValueError("Boundary must be placed at limit of the grid, indices 0 or -1")
                quit()
            else:
                self.slicing = slicing
                self.field = "E" if point_line == 0 else "H"
                self.orient()
        elif self.grid.dim == 3:
            types_list = [type(s) for s in slicing]
            if types_list.count(int) != 1:
                raise AssertionError("In three dimensions boundary is defined by a surface")
                quit()
            if not slicing[types_list.index(int)] in [0,-1]:
                raise ValueError("Boundary must be placed at limit of the grid, indices 0 or -1")
                quit()
            else:
                self.slicing = slicing
                self.field = "E" if slicing[types_list.index(int)] == 0 else "H"
                self.orient()
    
    def orient(self):
        is_flat_axis = [n ==1 for n in [self.grid.Nx, self.grid.Ny, self.grid.Nz]]
        orient_list = ["x","y","z"]
        if self.grid.dim == 1:
            self.orientation = orient_list[is_flat_axis.index(False)]
        elif self.grid.dim == 2:
            self.orientation = orient_list[[type(s) is int and not is_flat_axis[k] for k,s in enumerate(self.slicing)].index(True)]
        elif self.grid.dim == 3:
            self.orientation = orient_list[[type(s) for s in self.slicing].index(int)]
    
    def other_coefs(self):
        pass


class Periodic(Boundary):
    def update_H(self):
        if self.orientation == "x":
            self.grid.H[self.slicing][...,1] += self.grid.chye[self.slicing] * squeeze(self.grid.E[0,:,:,2] - self.grid.E[-1,:,:,2]) #Hy
            self.grid.H[self.slicing][...,2] -= self.grid.chze[self.slicing] * squeeze(self.grid.E[0,:,:,1] - self.grid.E[-1,:,:,1]) #Hz
        elif self.orientation == "y":
            self.grid.H[self.slicing][...,0] -= self.grid.chxe[self.slicing] * squeeze(self.grid.E[:,0,:,2] - self.grid.E[:,-1,:,2]) #Hx
            self.grid.H[self.slicing][...,2] += self.grid.chze[self.slicing] * squeeze(self.grid.E[:,0,:,0] - self.grid.E[:,-1,:,0]) #Hz
        elif self.orientation == "z":
            self.grid.H[self.slicing][...,0] += self.grid.chxe[self.slicing] * squeeze(self.grid.E[:,:,0,1] - self.grid.E[:,:,-1,1]) #Hx
            self.grid.H[self.slicing][...,1] -= self.grid.chye[self.slicing] * squeeze(self.grid.E[:,:,0,0] - self.grid.E[:,:,-1,0]) #Hy
    
    def update_E(self):
        if self.orientation == "x":
            self.grid.E[self.slicing][...,1] -= self.grid.ceyh[self.slicing] * squeeze(self.grid.H[0,:,:,2] - self.grid.H[-1,:,:,2]) #Ey
            self.grid.E[self.slicing][...,2] += self.grid.cezh[self.slicing] * squeeze(self.grid.H[0,:,:,1] - self.grid.H[-1,:,:,1]) #Ez
        elif self.orientation == "y":
            self.grid.E[self.slicing][...,0] += self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,0,:,2] - self.grid.H[:,-1,:,2]) #Ex
            self.grid.E[self.slicing][...,2] -= self.grid.cezh[self.slicing] * squeeze(self.grid.H[:,0,:,0] - self.grid.H[:,-1,:,0]) #Ez
        elif self.orientation == "z":
            self.grid.E[self.slicing][...,0] -= self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,:,0,1] - self.grid.H[:,:,-1,1]) #Ex
            self.grid.E[self.slicing][...,1] += self.grid.ceyh[self.slicing] * squeeze(self.grid.H[:,:,0,0] - self.grid.H[:,:,-1,0]) #Ey


class Mirror(Boundary):
    def update_H(self):
        if self.orientation == "x":
            self.grid.H[self.slicing][...,1] -= self.grid.chye[self.slicing] * squeeze(self.grid.E[-1,:,:,2]) #Hy
            self.grid.H[self.slicing][...,2] += self.grid.chze[self.slicing] * squeeze(self.grid.E[-1,:,:,1]) #Hz
        elif self.orientation == "y":
            self.grid.H[self.slicing][...,0] += self.grid.chxe[self.slicing] * squeeze(self.grid.E[:,-1,:,2]) #Hx
            self.grid.H[self.slicing][...,2] -= self.grid.chze[self.slicing] * squeeze(self.grid.E[:,-1,:,0]) #Hz
        elif self.orientation == "z":
            self.grid.H[self.slicing][...,0] -= self.grid.chxe[self.slicing] * squeeze(self.grid.E[:,:,-1,1]) #Hx
            self.grid.H[self.slicing][...,1] += self.grid.chye[self.slicing] * squeeze(self.grid.E[:,:,-1,0]) #Hy
    
    def update_E(self):
        if self.orientation == "x":
            self.grid.E[self.slicing][...,1] -= self.grid.ceyh[self.slicing] * squeeze(self.grid.H[0,:,:,2]) #Ey
            self.grid.E[self.slicing][...,2] += self.grid.cezh[self.slicing] * squeeze(self.grid.H[0,:,:,1]) #Ez
        elif self.orientation == "y":
            self.grid.E[self.slicing][...,0] += self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,0,:,2]) #Ex
            self.grid.E[self.slicing][...,2] -= self.grid.cezh[self.slicing] * squeeze(self.grid.H[:,0,:,0]) #Ez
        elif self.orientation == "z":
            self.grid.E[self.slicing][...,0] -= self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,:,0,1]) #Ex
            self.grid.E[self.slicing][...,1] += self.grid.ceyh[self.slicing] * squeeze(self.grid.H[:,:,0,0]) #Ey


class Bloch(Boundary):
    def __init__(self, bloch_vector):
        self.bloch_vector = bloch_vector)
        super().__init__()
    
    def other_coefs(self):
        #Check bloch vector is on grid
        is_flat_axis = [n==1 for n in [self.grid.Nx, self.grid.Ny, self.grid.Nz]]
        if any([self.bloch_vector[i]!=0 and not is_flat_axis[i] for i in range(3)]):
            raise AssertionError("Bloch vector out of dimensions' grid")
            quit()


class Absorbing(Boundary):
    def __init__(self):
        super().__init__()
        self.previous_H = None
        self.previous_E = None
        
        self.previous_extxHy = 0
        self.previous_extxHz = 0
        self.previous_extyHx = 0
        self.previous_extyHz = 0
        self.previous_extzHx = 0
        self.previous_extzHy = 0
        self.previous_extxEy = 0
        self.previous_extxEz = 0
        self.previous_extyEx = 0
        self.previous_extyEz = 0
        self.previous_extzEx = 0
        self.previous_extzEy = 0
        
        self.adv_coeff = None
    
    def other_coefs(self):
        reduced_courant = self.grid.courant_number/self.grid.refractive_index
        self.adv_coeff = (reduced_courant-1)/(reduced_courant+1)
        self.previous_H = copy(self.grid.H)
        self.previous_E = copy(self.grid.E)
    
    def update_H(self):
        if self.orientation == "x":
            extxEz = self.previous_E[-1,:,:,2] + self.adv_coeff[-1,:,:,2] * (self.grid.E[-1,:,:,2] - self.previous_extxEz)
            self.grid.H[self.slicing][...,1] += self.grid.chye[self.slicing] * squeeze(extxEz - self.grid.E[-1,:,:,2]) #Hy
            self.previous_extxEz = extxEz
            extxEy = self.previous_E[-1,:,:,1] + self.adv_coeff[-1,:,:,1] * (self.grid.E[-1,:,:,1] - self.previous_extxEy)
            self.grid.H[self.slicing][...,2] -= self.grid.chze[self.slicing] * squeeze(extxEy - self.grid.E[-1,:,:,1]) #Hz
            self.previous_extxEy = extxEy
        if self.orientation == "y":
            extyEz = self.previous_E[:,-1,:,2] + self.adv_coeff[:,-1,:,2] * (self.grid.E[:,-1,:,2] - self.previous_extyEz)
            self.grid.H[self.slicing][...,0] -= self.grid.chxe[self.slicing] * squeeze(extyEz - self.grid.E[:,-1,:,2]) #Hx
            self.previous_extyEz = extyEz
            extyEx = self.previous_E[:,-1,:,0] + self.adv_coeff[:,-1,:,0] * (self.grid.E[:,-1,:,0] - self.previous_extyEx)
            self.grid.H[self.slicing][...,2] += self.grid.chze[self.slicing] * squeeze(extyEx - self.grid.E[:,-1,:,0]) #Hz
            self.previous_extyEx = extyEx
        if self.orientation == "z":
            extzEy = self.previous_E[:,:,-1,1] + self.adv_coeff[:,:,-1,1] * (self.grid.E[:,:,-1,1] - self.previous_extzEy)
            self.grid.H[self.slicing][...,0] += self.grid.chxe[self.slicing] * squeeze(extzEy - self.grid.E[:,:,-1,1]) #Hx
            self.previous_extzEy = extzEy
            extzEx = self.previous_E[:,:,-1,0] + self.adv_coeff[:,:,-1,0] * (self.grid.E[:,:,-1,0] - self.previous_extzEx)
            self.grid.H[self.slicing][...,1] -= self.grid.chye[self.slicing] * squeeze(extzEx - self.grid.E[:,:,-1,0]) #Hy
            self.previous_extzEx = extzEx
        self.previous_E = copy(self.grid.E)
    
    def update_E(self):
        if self.orientation == "x":
            extxHz = self.previous_H[0,:,:,2] + self.adv_coeff[0,:,:,2] * (self.grid.H[0,:,:,2] - self.previous_extxHz)
            self.grid.E[self.slicing][...,1] -= self.grid.ceyh[self.slicing] * squeeze(self.grid.H[0,:,:,2] - extxHz) #Ey
            self.previous_extxHz = extxHz
            extxHy = self.previous_H[0,:,:,1] + self.adv_coeff[0,:,:,1] * (self.grid.H[0,:,:,1] - self.previous_extxHy)
            self.grid.E[self.slicing][...,2] += self.grid.cezh[self.slicing] * squeeze(self.grid.H[0,:,:,1] - extxHy) #Ez
            self.previous_extxHy = extxHy
        if self.orientation == "y":
            extyHz = self.previous_H[:,0,:,2] + self.adv_coeff[:,0,:,2] * (self.grid.H[:,0,:,2] - self.previous_extyHz)
            self.grid.E[self.slicing][...,0] += self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,0,:,2] - extyHz) #Ex
            self.previous_extyHz = extyHz
            extyHx = self.previous_H[:,0,:,0] + self.adv_coeff[:,0,:,0] * (self.grid.H[:,0,:,0] - self.previous_extyHx)
            self.grid.E[self.slicing][...,2] -= self.grid.cezh[self.slicing] * squeeze(self.grid.H[:,0,:,0] - extyHx) #Ez
            self.previous_extyHx = extyHx
        if self.orientation == "z":
            extzHy = self.previous_H[:,:,0,1] + self.adv_coeff[:,:,0,1] * (self.grid.H[:,:,0,1] - self.previous_extzHy)
            self.grid.E[self.slicing][...,0] -= self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,:,0,1] - extzHy) #Ex
            self.previous_extzHy = extzHy
            extzHx = self.previous_H[:,:,0,0] + self.adv_coeff[:,:,0,0] * (self.grid.H[:,:,0,0] - self.previous_extzHx)
            self.grid.E[self.slicing][...,1] += self.grid.ceyh[self.slicing] * squeeze(self.grid.H[:,:,0,0] - extzHx) #Ey
            self.previous_extzHx = extzHx
        self.previous_H = copy(self.grid.H)
