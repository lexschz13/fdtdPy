from pylab import *
from .vacuum_constants import *



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
    
    def update_boundH(self):
        pass
    
    def update_boundE(self):
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
        self.bloch_vector = bloch_vector
        super().__init__()
    
    def other_coefs(self):
        #Check bloch vector is on grid
        is_flat_axis = [n==1 for n in [self.grid.Nx, self.grid.Ny, self.grid.Nz]]
        if any([self.bloch_vector[i]!=0 and is_flat_axis[i] for i in range(3)]):
            raise AssertionError("Bloch vector out of dimensions' grid")
            quit()
        #Check if grid is using complex numbers
        if not self.grid.is_phasor:
            raise AssertionError("Bloch boundary condition only can be added if the simulation is using complex numbers")
            quit()
    
    def update_H(self):
        if self.orientation == "x":
            a = self.grid.spatial_step * self.grid.Nx
            c = exp(-1j * (a*self.bloch_vector[0])) #Correction bloch factor
            self.grid.H[self.slicing][...,1] += self.grid.chye[self.slicing] * squeeze(c*self.grid.E[0,:,:,2] - self.grid.E[-1,:,:,2]) #Hy
            self.grid.H[self.slicing][...,2] -= self.grid.chze[self.slicing] * squeeze(c*self.grid.E[0,:,:,1] - self.grid.E[-1,:,:,1]) #Hz
        elif self.orientation == "y":
            a = self.grid.spatial_step * self.grid.Ny
            c = exp(-1j * (a*self.bloch_vector[1])) #Correction bloch factor
            self.grid.H[self.slicing][...,0] -= self.grid.chxe[self.slicing] * squeeze(c*self.grid.E[:,0,:,2] - self.grid.E[:,-1,:,2]) #Hx
            self.grid.H[self.slicing][...,2] += self.grid.chze[self.slicing] * squeeze(c*self.grid.E[:,0,:,0] - self.grid.E[:,-1,:,0]) #Hz
        elif self.orientation == "z":
            a = self.grid.spatial_step * self.grid.Nz
            c = exp(-1j * (a*self.bloch_vector[2])) #Correction bloch factor
            self.grid.H[self.slicing][...,0] += self.grid.chxe[self.slicing] * squeeze(c*self.grid.E[:,:,0,1] - self.grid.E[:,:,-1,1]) #Hx
            self.grid.H[self.slicing][...,1] -= self.grid.chye[self.slicing] * squeeze(c*self.grid.E[:,:,0,0] - self.grid.E[:,:,-1,0]) #Hy
    
    def update_E(self):
        if self.orientation == "x":
            a = self.grid.spatial_step * self.grid.Nx
            c = exp(+1j * (a*self.bloch_vector[0])) #Correction bloch factor
            self.grid.E[self.slicing][...,1] -= self.grid.ceyh[self.slicing] * squeeze(self.grid.H[0,:,:,2] - c*self.grid.H[-1,:,:,2]) #Ey
            self.grid.E[self.slicing][...,2] += self.grid.cezh[self.slicing] * squeeze(self.grid.H[0,:,:,1] - c*self.grid.H[-1,:,:,1]) #Ez
        elif self.orientation == "y":
            a = self.grid.spatial_step * self.grid.Ny
            c = exp(+1j * (a*self.bloch_vector[1])) #Correction bloch factor
            self.grid.E[self.slicing][...,0] += self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,0,:,2] - c*self.grid.H[:,-1,:,2]) #Ex
            self.grid.E[self.slicing][...,2] -= self.grid.cezh[self.slicing] * squeeze(self.grid.H[:,0,:,0] - c*self.grid.H[:,-1,:,0]) #Ez
        elif self.orientation == "z":
            a = self.grid.spatial_step * self.grid.Nz
            c = exp(+1j * (a*self.bloch_vector[2])) #Correction bloch factor
            self.grid.E[self.slicing][...,0] -= self.grid.cexh[self.slicing] * squeeze(self.grid.H[:,:,0,1] - c*self.grid.H[:,:,-1,1]) #Ex
            self.grid.E[self.slicing][...,1] += self.grid.ceyh[self.slicing] * squeeze(self.grid.H[:,:,0,0] - c*self.grid.H[:,:,-1,0]) #Ey


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


class PML(Boundary):
    def __init__(self, refractive_index=1, kappa=1, alpha=1, sigma=1, thickness=10, growing_step=1e-7):
        super().__init__()
        self.previous_E = None #For own abs bounds
        self.previous_H = None
        self.adv_coeff = None
        
        self.kappa = kappa
        self.alpha = alpha
        self.sigma = sigma
        
        self.growing_step = growing_step
        
        if refractive_index < 1:
            raise ValueError("Refractive index has to be 1 or more")
            quit()
        else:
            self.refractive_index = refractive_index
        
        if thickness <= 1:
            raise ValueError("Thickess has to be at least 2")
            quit()
        else:
            self.thickness = thickness
        
        self.Nx, self.Ny, self.Nz = [None]*3
        
        self.b = None
        self.C = None
        
        self.E = None
        self.H = None
        
        self.psiH = None
        self.psiE = None
        
        self.chxe, self.chye, self.chze, self.cexh, self.ceyh, self.cezh = [None]*6
        self.adv_coeff = None
        
        self.prev_back_extxEy = None
        self.prev_back_extxEz = None
        self.prev_back_extyEx = None
        self.prev_back_extyEz = None
        self.prev_back_extzEx = None
        self.prev_back_extzEy = None
        
        self.prev_trans_extxEy = None
        self.prev_trans_extxEz = None
        self.prev_trans_extyEx = None
        self.prev_trans_extyEz = None
        self.prev_trans_extzEx = None
        self.prev_trans_extzEy = None
        
        self.prev_back_extxHy = None
        self.prev_back_extxHz = None
        self.prev_back_extyHx = None
        self.prev_back_extyHz = None
        self.prev_back_extzHx = None
        self.prev_back_extzHy = None
        
        self.prev_trans_extxHy = None
        self.prev_trans_extxHz = None
        self.prev_trans_extyHx = None
        self.prev_trans_extyHz = None
        self.prev_trans_extzHx = None
        self.prev_trans_extzHy = None
        
    def set_coefs(self):
        #mu1/eps1 = mu2/eps2
        #nPML = sqrt(eps2*mu2)*c = eps2 * sqrt(mu1/eps1) * c = eps2/eps1 * ng  => eps2 = eps1 * nPML/ng  and  mu2 = mu1 * nPML/ng
        
        #Unsqueezing
        if self.grid.dim == 1:
            flat_sh = (1, 1, 1, 3)
        elif self.grid.dim == 2:
            if self.orientation == "x":
                flat_sh = (1, 1, self.Nz, 3) if self.grid.Ny == 1 else (1, self.Ny, 1, 3)
            elif self.orientation == "y":
                flat_sh = (self.Nx, 1, 1, 3) if self.grid.Nz == 1 else (1, 1, self.Nz, 3)
            elif self.orientation == "z":
                flat_sh = (1, self.Ny, 1, 3) if self.grid.Nx == 1 else (self.Nx, 1, 1, 3)
        elif self.grid.dim == 3:
            if self.orientation == "x":
                flat_sh = (1, self.Ny, self.Nz, 3)
            elif self.orientation == "y":
                flat_sh = (self.Nx, 1, self.Nz, 3)
            elif self.orientation == "z":
                flat_sh = (self.Nx, self.Ny, 1, 3)
        
        ng = self.grid.refractive_index[self.slicing].reshape(flat_sh)
        epsg = self.grid.permittivity[self.slicing].reshape(flat_sh)
        mug = self.grid.permeability[self.slicing].reshape(flat_sh)
        self.refractive_index *= ones((self.Nx, self.Ny, self.Nz, 3))
        epsPML = epsg * self.refractive_index/ng
        muPML = mug * self.refractive_index/ng
        
        self.chxe = self.grid.courant_number / muPML[...,0] / VACUUM_IMPEDANCE
        self.chye = self.grid.courant_number / muPML[...,1] / VACUUM_IMPEDANCE
        self.chze = self.grid.courant_number / muPML[...,2] / VACUUM_IMPEDANCE
        self.cexh = self.grid.courant_number / epsPML[...,0] * VACUUM_IMPEDANCE
        self.ceyh = self.grid.courant_number / epsPML[...,1] * VACUUM_IMPEDANCE
        self.cezh = self.grid.courant_number / epsPML[...,2] * VACUUM_IMPEDANCE
        
        reduced_courant = self.grid.courant_number/self.refractive_index
        self.adv_coeff = (reduced_courant-1)/(reduced_courant+1)
    
    def back_extxEy(self, update_adv=False):
        #Background boundary-x
        extxEy = self.previous_E[-1,:,:,1] + self.adv_coeff[-1,:,:,1] * (self.E[-1,:,:,1] - self.prev_back_extxEy)
        if update_adv: self.prev_back_extxEy = extxEy
        return extxEy
    
    def back_extxEz(self, update_adv=False):
        #Background boundary-x
        extxEz = self.previous_E[-1,:,:,2] + self.adv_coeff[-1,:,:,2] * (self.E[-1,:,:,2] - self.prev_back_extxEz)
        if update_adv: self.prev_back_extxEz = extxEz
        return extxEz
    
    def back_extyEx(self, update_adv=False):
        #Background boundary-y
        extyEx = self.previous_E[:,-1,:,0] + self.adv_coeff[:,-1,:,0] * (self.E[:,-1,:,0] - self.prev_back_extyEx)
        if update_adv: self.prev_back_extyEx = extyEx
        return extyEx
    
    def back_extyEz(self, update_adv=False):
        #Background boundary-y
        extyEz = self.previous_E[:,-1,:,2] + self.adv_coeff[:,-1,:,2] * (self.E[:,-1,:,2] - self.prev_back_extyEz)
        if update_adv: self.prev_back_extyEz = extyEz
        return extyEz
    
    def back_extzEx(self, update_adv=False):
        #Background boundary-z
        extzEx = self.previous_E[:,:,-1,0] + self.adv_coeff[:,:,-1,0] * (self.E[:,:,-1,0] - self.prev_back_extzEx)
        if update_adv: self.prev_back_extzEx = extzEx
        return extzEx
    
    def back_extzEy(self, update_adv=False):
        #Background boundary-z
        extzEy = self.previous_E[:,:,-1,1] + self.adv_coeff[:,:,-1,1] * (self.E[:,:,-1,1] - self.prev_back_extzEy)
        if update_adv: self.prev_back_extzEy = extzEy
        return extzEy
    
    def back_extxHy(self, update_adv=False):
        #Background boundary-x
        extxHy = self.previous_H[0,:,:,1] + self.adv_coeff[0,:,:,1] * (self.H[0,:,:,1] - self.prev_back_extxHy)
        if update_adv: self.prev_back_extxHy = extxHy
        return extxHy
    
    def back_extxHz(self, update_adv=False):
        #Background boundary-x
        extxHz = self.previous_H[0,:,:,2] + self.adv_coeff[0,:,:,2] * (self.H[0,:,:,2] - self.prev_back_extxHz)
        if update_adv: self.prev_back_extxHz = extxHz
        return extxHz
    
    def back_extyHx(self, update_adv=False):
        #Background boundary-y
        extyHx = self.previous_H[:,0,:,0] + self.adv_coeff[:,0,:,0] * (self.H[:,0,:,0] - self.prev_back_extyHx)
        if update_adv: self.prev_back_extyHx = extyHx
        return extyHx
    
    def back_extyHz(self, update_adv=False):
        #Background boundary-y
        extyHz = self.previous_H[:,0,:,2] + self.adv_coeff[:,0,:,2] * (self.H[:,0,:,2] - self.prev_back_extyHz)
        if update_adv: self.prev_back_extyHz = extyHz
        return extyHz
    
    def back_extzHx(self, update_adv=False):
        #Background boundary-z
        extzHx = self.previous_H[:,:,0,0] + self.adv_coeff[:,:,0,0] * (self.H[:,:,0,0] - self.prev_back_extzHx)
        if update_adv: self.prev_back_extzHx = extzHx
        return extzHx
    
    def back_extzHy(self, update_adv=False):
        #Background boundary-z
        extzHy = self.previous_H[:,:,0,1] + self.adv_coeff[:,:,0,1] * (self.H[:,:,0,1] - self.prev_back_extzHy)
        if update_adv: self.prev_back_extzHy = extzHy
        return extzHy
    
    def trans_extxEy(self, update_adv=False):
        #Transversal boundary-x
        if isinstance(self.grid.boundaries["xH"], Periodic):
            return self.E[0,:,:,1]
        elif isinstance(self.grid.boundaries["xH"], Bloch):
            a = self.grid.spatial_step * self.Nx
            c = exp(-1j * (a*self.grid.boundaries["xH"].bloch_vector[0])) #Correction bloch factor
            return c * self.E[0,:,:,1]
        elif isinstance(self.grid.boundaries["xH"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["xH"], (Absorbing, PML)):
            extxEy = self.previous_E[-1,:,:,1] + self.adv_coeff[-1,:,:,1] * (self.E[-1,:,:,1] - self.prev_trans_extxEy)
            if update_adv: self.prev_trans_extxEy = extxEy
            return extxEy
    
    def trans_extxEz(self, update_adv=False):
        #Transversal boundary-x
        if isinstance(self.grid.boundaries["xH"], Periodic):
            return self.E[0,:,:,2]
        elif isinstance(self.grid.boundaries["xH"], Bloch):
            a = self.grid.spatial_step * self.Nx
            c = exp(-1j * (a*self.grid.boundaries["xH"].bloch_vector[0])) #Correction bloch factor
            return c * self.E[0,:,:,2]
        elif isinstance(self.grid.boundaries["xH"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["xH"], (Absorbing, PML)):
            extxEz = self.previous_E[-1,:,:,2] + self.adv_coeff[-1,:,:,2] * (self.E[-1,:,:,2] - self.prev_trans_extxEz)
            if update_adv: self.prev_trans_extxEz = extxEz
            return extxEz
    
    def trans_extyEx(self, update_adv=False):
        #Transversal boundary-y
        if isinstance(self.grid.boundaries["yH"], Periodic):
            return self.E[:,0,:,0]
        elif isinstance(self.grid.boundaries["yH"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(-1j * (a*self.grid.boundaries["yH"].bloch_vector[1])) #Correction bloch factor
            return c * self.E[:,0,:,0]
        elif isinstance(self.grid.boundaries["yH"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["yH"], (Absorbing, PML)):
            extyEx = self.previous_E[:,-1,:,0] + self.adv_coeff[:,-1,:,0] * (self.E[:,-1,:,0] - self.prev_trans_extyEx)
            if update_adv: self.prev_trans_extyEx = extyEx
            return extyEx
    
    def trans_extyEz(self, update_adv=False):
        #Transversal boundary-y
        if isinstance(self.grid.boundaries["yH"], Periodic):
            return self.E[:,0,:,2]
        elif isinstance(self.grid.boundaries["yH"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(-1j * (a*self.grid.boundaries["yH"].bloch_vector[1])) #Correction bloch factor
            return c * self.E[:,0,:,2]
        elif isinstance(self.grid.boundaries["yH"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["yH"], (Absorbing, PML)):
            extyEz = self.previous_E[:,-1,:,2] + self.adv_coeff[:,-1,:,2] * (self.E[:,-1,:,2] - self.prev_trans_extyEz)
            if update_adv: self.prev_trans_extyEz = extyEz
            return extyEz
    
    def trans_extzEx(self, update_adv=False):
        #Transversal boundary-z
        if isinstance(self.grid.boundaries["zH"], Periodic):
            return self.E[:,:,0,0]
        elif isinstance(self.grid.boundaries["zH"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(-1j * (a*self.grid.boundaries["zH"].bloch_vector[2])) #Correction bloch factor
            return c * self.E[:,:,0,0]
        elif isinstance(self.grid.boundaries["zH"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["zH"], (Absorbing, PML)):
            extzEx = self.previous_E[:,:,-1,0] + self.adv_coeff[:,:,-1,0] * (self.E[:,:,-1,0] - self.prev_trans_extzEx)
            if update_adv: self.prev_trans_extzEx = extzEx
            return extzEx
    
    def trans_extzEy(self, update_adv=False):
        #Transversal boundary-z
        if isinstance(self.grid.boundaries["zH"], Periodic):
            return self.E[:,:,0,1]
        elif isinstance(self.grid.boundaries["zH"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(-1j * (a*self.grid.boundaries["zH"].bloch_vector[2])) #Correction bloch factor
            return c * self.E[:,:,0,1]
        elif isinstance(self.grid.boundaries["zH"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["zH"], (Absorbing, PML)):
            extzEy = self.previous_E[:,:,-1,1] + self.adv_coeff[:,:,-1,1] * (self.E[:,:,-1,1] - self.prev_trans_extzEy)
            if update_adv: self.prev_trans_extzEy = extzEy
            return extzEy
    
    def trans_extxHy(self, update_adv=False):
        #Transversal boundary-x
        if isinstance(self.grid.boundaries["xE"], Periodic):
            return self.H[-1,:,:,1]
        elif isinstance(self.grid.boundaries["xE"], Bloch):
            a = self.grid.spatial_step * self.Nx
            c = exp(+1j * (a*self.grid.boundaries["xE"].bloch_vector[0])) #Correction bloch factor
            return c * self.H[-1,:,:,1]
        elif isinstance(self.grid.boundaries["xE"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["xE"], (Absorbing, PML)):
            extxHy = self.previous_H[0,:,:,1] + self.adv_coeff[0,:,:,1] * (self.H[0,:,:,1] - self.prev_trans_extxHy)
            if update_adv: self.prev_trans_extxHy = extxHy
            return extxHy
    
    def trans_extxHz(self, update_adv=False):
        #Transversal boundary-x
        if isinstance(self.grid.boundaries["xE"], Periodic):
            return self.H[-1,:,:,2]
        elif isinstance(self.grid.boundaries["xE"], Bloch):
            a = self.grid.spatial_step * self.Nx
            c = exp(+1j * (a*self.grid.boundaries["xE"].bloch_vector[0])) #Correction bloch factor
            return c * self.H[-1,:,:,2]
        elif isinstance(self.grid.boundaries["xE"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["xE"], (Absorbing, PML)):
            extxHz = self.previous_H[0,:,:,2] + self.adv_coeff[0,:,:,2] * (self.H[0,:,:,2] - self.prev_trans_extxHz)
            if update_adv: self.prev_trans_extxHz = extxHz
            return extxHz
    
    def trans_extyHx(self, update_adv=False):
        #Transversal boundary-y
        if isinstance(self.grid.boundaries["yE"], Periodic):
            return self.H[:,-1,:,0]
        elif isinstance(self.grid.boundaries["yE"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(+1j * (a*self.grid.boundaries["yE"].bloch_vector[1])) #Correction bloch factor
            return c * self.H[:,-1,:,0]
        elif isinstance(self.grid.boundaries["yE"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["yE"], (Absorbing, PML)):
            extyHx = self.previous_H[:,0,:,0] + self.adv_coeff[:,0,:,0] * (self.H[:,0,:,0] - self.prev_trans_extyHx)
            if update_adv: self.prev_trans_extyHx = extyHx
            return extyHx
    
    def trans_extyHz(self, update_adv=False):
        #Transversal boundary-y
        if isinstance(self.grid.boundaries["yE"], Periodic):
            return self.H[:,-1,:,2]
        elif isinstance(self.grid.boundaries["yE"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(+1j * (a*self.grid.boundaries["yE"].bloch_vector[1])) #Correction bloch factor
            return c * self.H[:,-1,:,2]
        elif isinstance(self.grid.boundaries["yE"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["yE"], (Absorbing, PML)):
            extyHz = self.previous_H[:,0,:,2] + self.adv_coeff[:,0,:,2] * (self.H[:,0,:,2] - self.prev_trans_extyHz)
            if update_adv: self.prev_trans_extyHz = extyHz
            return extyHz
    
    def trans_extzHx(self, update_adv=False):
        #Transversal boundary-z
        if isinstance(self.grid.boundaries["zE"], Periodic):
            return self.H[:,:,-1,0]
        elif isinstance(self.grid.boundaries["zE"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(+1j * (a*self.grid.boundaries["zE"].bloch_vector[2])) #Correction bloch factor
            return c * self.H[:,:,-1,0]
        elif isinstance(self.grid.boundaries["zE"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["zE"], (Absorbing, PML)):
            extzHx = self.previous_H[:,:,0,0] + self.adv_coeff[:,:,0,0] * (self.H[:,:,0,0] - self.prev_trans_extzHx)
            if update_adv: self.prev_trans_extzHx = extzHx
            return extzHx
    
    def trans_extzHy(self, update_adv=False):
        #Transversal boundary-z
        if isinstance(self.grid.boundaries["zE"], Periodic):
            return self.H[:,:,-1,1]
        elif isinstance(self.grid.boundaries["zE"], Bloch):
            a = self.grid.spatial_step * self.Ny
            c = exp(+1j * (a*self.grid.boundaries["zE"].bloch_vector[2])) #Correction bloch factor
            return c * self.H[:,:,-1,1]
        elif isinstance(self.grid.boundaries["zE"], Mirror):
            return 0
        elif isinstance(self.grid.boundaries["zE"], (Absorbing, PML)):
            extzHy = self.previous_H[:,:,0,1] + self.adv_coeff[:,:,0,1] * (self.H[:,:,0,1] - self.prev_trans_extzHy)
            if update_adv: self.prev_trans_extzHy = extzHy
            return extzHy
    
    def poly(self):
        #Defines a cubic growth
        if self.field == "H":
            x = arange(0,self.thickness,1)*self.grid.spatial_step/self.growing_step
        elif self.field == "E":
            x = arange(self.thickness-1,-1,-1)*self.grid.spatial_step/self.growing_step
        return x**3
    
    def other_coefs(self):
        if self.orientation == "x":
            self.kappa = array([self.kappa,1,1])[None,None,None,:]
            self.sigma = array([self.sigma,0,0])[None,None,None,:] * self.poly()[:,None,None,None]
            self.Nx, self.Ny, self.Nz = self.thickness, self.grid.Ny, self.grid.Nz
            self.prev_back_extxEy = zeros((self.Ny,self.Nz))
            self.prev_back_extxEz = zeros((self.Ny,self.Nz))
            self.prev_trans_extyEx = zeros((self.Nx,self.Nz))
            self.prev_trans_extyEz = zeros((self.Nx,self.Nz))
            self.prev_trans_extzEx = zeros((self.Nx,self.Ny))
            self.prev_trans_extzEy = zeros((self.Nx,self.Ny))
            self.prev_back_extxHy = zeros((self.Ny,self.Nz))
            self.prev_back_extxHz = zeros((self.Ny,self.Nz))
            self.prev_trans_extyHx = zeros((self.Nx,self.Nz))
            self.prev_trans_extyHz = zeros((self.Nx,self.Nz))
            self.prev_trans_extzHx = zeros((self.Nx,self.Ny))
            self.prev_trans_extzHy = zeros((self.Nx,self.Ny))
        elif self.orientation == "y":
            self.kappa = array([1,self.kappa,1])[None,None,None,:]
            self.sigma = array([0,self.sigma,0])[None,None,None,:] * self.poly()[None,:,None,None]
            self.Nx, self.Ny, self.Nz = self.grid.Nx, self.thickness, self.grid.Nz
            self.prev_back_extyEx = zeros((self.Nx,self.Nz))
            self.prev_back_extyEz = zeros((self.Nx,self.Nz))
            self.prev_trans_extxEy = zeros((self.Ny,self.Nz))
            self.prev_trans_extxEz = zeros((self.Ny,self.Nz))
            self.prev_trans_extzEx = zeros((self.Nx,self.Ny))
            self.prev_trans_extzEy = zeros((self.Nx,self.Ny))
            self.prev_back_extyHx = zeros((self.Nx,self.Nz))
            self.prev_back_extyHz = zeros((self.Nx,self.Nz))
            self.prev_trans_extxHy = zeros((self.Ny,self.Nz))
            self.prev_trans_extxHz = zeros((self.Ny,self.Nz))
            self.prev_trans_extzHx = zeros((self.Nx,self.Ny))
            self.prev_trans_extzHy = zeros((self.Nx,self.Ny))
        elif self.orientation == "z":
            self.kappa = array([1,1,self.kappa])[None,None,None,:]
            self.sigma = array([0,0,self.sigma])[None,None,None,:] * self.poly()[None,None,:,None]
            self.Nx, self.Ny, self.Nz = self.grid.Nx, self.grid.Ny, self.thickness
            self.prev_back_extzEx = zeros((self.Nx,self.Ny))
            self.prev_back_extzEy = zeros((self.Nx,self.Ny))
            self.prev_trans_extxEy = zeros((self.Ny,self.Nz))
            self.prev_trans_extxEz = zeros((self.Ny,self.Nz))
            self.prev_trans_extyEx = zeros((self.Nx,self.Nz))
            self.prev_trans_extyEz = zeros((self.Nx,self.Nz))
            self.prev_back_extzHx = zeros((self.Nx,self.Ny))
            self.prev_back_extzHy = zeros((self.Nx,self.Ny))
            self.prev_trans_extxHy = zeros((self.Ny,self.Nz))
            self.prev_trans_extxHz = zeros((self.Ny,self.Nz))
            self.prev_trans_extyHx = zeros((self.Nx,self.Nz))
            self.prev_trans_extyHz = zeros((self.Nx,self.Nz))
        self.b = exp(-self.grid.time_pace/VACUUM_PERMITTIVITY * (self.alpha+self.sigma/self.kappa)) * ones((self.Nx,self.Ny,self.Nz,3))
        self.C = self.sigma/self.kappa / (self.sigma + self.kappa*self.alpha) * (self.b-1)
        F_TYPE = complex128 if self.grid.is_phasor else float64
        self.H = zeros((self.Nx, self.Ny, self.Nz, 3), dtype=F_TYPE)
        self.E = zeros((self.Nx, self.Ny, self.Nz, 3), dtype=F_TYPE)
        self.previous_H = copy(self.H)
        self.previous_E = copy(self.E)
        self.psiH = zeros((self.Nx, self.Ny, self.Nz, 3, 3), dtype=F_TYPE) #Shape of layer grid, derivative direction, field component
        self.psiE = zeros((self.Nx, self.Ny, self.Nz, 3, 3), dtype=F_TYPE)
        self.set_coefs()
    
    def update_psiH(self):
        self.psiH *= self.b[...,None]
        if self.Nx > 1:
            self.psiH[:-1,:,:,0,1] += self.C[:-1,:,:,0] * (self.E[1:,:,:,2] - self.E[:-1,:,:,2])
            self.psiH[:-1,:,:,0,2] += self.C[:-1,:,:,0] * (self.E[1:,:,:,1] - self.E[:-1,:,:,1])
            #Boundaries
            if self.orientation == "x": #Normal boundaries
                if self.field == "E": #At low extreme of grid, -1 idx is joined to grid
                    self.psiH[-1,:,:,0,1] += self.C[-1,:,:,0] * (self.grid.E[0,:,:,2] - self.E[-1,:,:,2])
                    self.psiH[-1,:,:,0,2] += self.C[-1,:,:,0] * (self.grid.E[0,:,:,1] - self.E[-1,:,:,1])
                elif self.field == "H": #At high extreme of grid
                    self.psiH[-1,:,:,0,1] += self.C[-1,:,:,0] * (self.back_extxEz(False) - self.E[-1,:,:,2])
                    self.psiH[-1,:,:,0,2] += self.C[-1,:,:,0] * (self.back_extxEy(False) - self.E[-1,:,:,1])
            else: #Transversal boundaries
                self.psiH[-1,:,:,0,1] += self.C[-1,:,:,0] * (self.trans_extxEz(False) - self.E[-1,:,:,2])
                self.psiH[-1,:,:,0,2] += self.C[-1,:,:,0] * (self.trans_extxEy(False) - self.E[-1,:,:,1])
        if self.Ny > 1:
            self.psiH[:,:-1,:,1,0] += self.C[:,:-1,:,1] * (self.E[:,1:,:,2] - self.E[:,:-1,:,2])
            self.psiH[:,:-1,:,1,2] += self.C[:,:-1,:,1] * (self.E[:,1:,:,0] - self.E[:,:-1,:,0])
            #Boundaries
            if self.orientation == "y": #Normal boundaries
                if self.field == "E": #At low extreme of grid, -1 idx is joined to grid
                    self.psiH[:,-1,:,1,0] += self.C[:,-1,:,1] * (self.grid.E[:,0,:,2] - self.E[:,-1,:,2])
                    self.psiH[:,-1,:,1,2] += self.C[:,-1,:,1] * (self.grid.E[:,0,:,0] - self.E[:,-1,:,0])
                elif self.field == "H": #At high extreme of grid
                    self.psiH[:,-1,:,1,0] += self.C[:,-1,:,1] * (self.back_extyEz(False) - self.E[:,-1,:,2])
                    self.psiH[:,-1,:,1,2] += self.C[:,-1,:,1] * (self.back_extyEx(False) - self.E[:,-1,:,0])
            else: #Transveral boundaries
                self.psiH[:,-1,:,1,0] += self.C[:,-1,:,1] * (self.trans_extyEz(False) - self.E[:,-1,:,2])
                self.psiH[:,-1,:,1,2] += self.C[:,-1,:,1] * (self.trans_extyEx(False) - self.E[:,-1,:,0])
        if self.Nz > 1:
            self.psiH[:,:,:-1,2,0] += self.C[:,:,:-1,2] * (self.E[:,:,1:,1] - self.E[:,:,:-1,1])
            self.psiH[:,:,:-1,2,1] += self.C[:,:,:-1,2] * (self.E[:,:,1:,0] - self.E[:,:,:-1,0])
            #Boundaries
            if self.orientation == "z": #Normal boundaries
                if self.field == "E": #At low extreme of grid, -1 idx is joined to grid
                    self.psiH[:,:,-1,2,0] += self.C[:,:,-1,2] * (self.grid.E[:,:,0,1] - self.E[:,:,-1,1])
                    self.psiH[:,:,-1,2,1] += self.C[:,:,-1,2] * (self.grid.E[:,:,0,0] - self.E[:,:,-1,0])
                elif self.field == "H": #At high extreme of grid
                    self.psiH[:,:,-1,2,0] += self.C[:,:,-1,2] * (self.back_extzEy(False) - self.E[:,:,-1,1])
                    self.psiH[:,:,-1,2,1] += self.C[:,:,-1,2] * (self.back_extzEx(False) - self.E[:,:,-1,0])
            else: #Transversal boundaries
                self.psiH[:,:,-1,2,0] += self.C[:,:,-1,2] * (self.trans_extzEy(False) - self.E[:,:,-1,1])
                self.psiH[:,:,-1,2,1] += self.C[:,:,-1,2] * (self.trans_extzEx(False) - self.E[:,:,-1,0])
    
    def update_boundH(self):
        self.update_psiH()
        #Update on x axis
        if self.Nx > 1:
            self.H[:-1,:,:,1] += self.chye[:-1,:,:] * ((self.E[1:,:,:,2] - self.E[:-1,:,:,2])/squeeze(self.kappa[...,0]) + self.psiH[:-1,:,:,0,1]) #Hy
            self.H[:-1,:,:,2] -= self.chze[:-1,:,:] * ((self.E[1:,:,:,1] - self.E[:-1,:,:,1])/squeeze(self.kappa[...,0]) + self.psiH[:-1,:,:,0,2]) #Hz
            #Boundaries
            if self.orientation == "x": #Normal boundaries
                if self.field == "E": #At low extreme of grid, -1 idx is joined to grid
                    self.H[-1,:,:,1] += self.chye[-1,:,:] * ((self.grid.E[0,:,:,2] - self.E[-1,:,:,2])/squeeze(self.kappa[...,0]) + self.psiH[-1,:,:,0,1]) #Hy
                    self.H[-1,:,:,2] -= self.chze[-1,:,:] * ((self.grid.E[0,:,:,1] - self.E[-1,:,:,1])/squeeze(self.kappa[...,0]) + self.psiH[-1,:,:,0,2]) #Hz
                elif self.field == "H": #At high extreme of grid
                    self.H[-1,:,:,1] += self.chye[-1,:,:] * ((self.back_extxEz(True) - self.E[-1,:,:,2])/squeeze(self.kappa[...,0]) + self.psiH[-1,:,:,0,1]) #Hy
                    self.H[-1,:,:,2] -= self.chze[-1,:,:] * ((self.back_extxEy(True) - self.E[-1,:,:,1])/squeeze(self.kappa[...,0]) + self.psiH[-1,:,:,0,2]) #Hz
            else: #Transversal boundaries
                self.H[-1,:,:,1] += self.chye[-1,:,:] * ((self.trans_extxEz(True) - self.E[-1,:,:,2])/squeeze(self.kappa[...,0]) + self.psiH[-1,:,:,0,1]) #Hy
                self.H[-1,:,:,2] -= self.chze[-1,:,:] * ((self.trans_extxEy(True) -self.E[-1,:,:,1])/squeeze(self.kappa[...,0]) + self.psiH[-1,:,:,0,2]) #Hz
        
        #Update on y axis
        if self.Ny > 1:
            self.H[:,:-1,:,0] -= self.chxe[:,:-1,:] * ((self.E[:,1:,:,2] - self.E[:,:-1,:,2])/squeeze(self.kappa[...,1]) + self.psiH[:,:-1,:,1,0]) #Hx
            self.H[:,:-1,:,2] += self.chze[:,:-1,:] * ((self.E[:,1:,:,0] - self.E[:,:-1,:,0])/squeeze(self.kappa[...,1]) + self.psiH[:,:-1,:,1,2]) #Hz
            #Boundaries
            if self.orientation == "y": #Normal boundaries
                if self.field == "E": #At low extreme of grid, -1 idx is joined to grid
                    self.H[:,-1,:,0] -= self.chxe[:,-1,:] * ((self.grid.E[:,0,:,2] - self.E[:,-1,:,2])/squeeze(self.kappa[...,1]) + self.psiH[:,-1,:,1,0]) #Hx
                    self.H[:,-1,:,2] += self.chze[:,-1,:] * ((self.grid.E[:,0,:,0] - self.E[:,-1,:,0])/squeeze(self.kappa[...,1]) + self.psiH[:,-1,:,1,2]) #Hz
                elif self.field == "H": #At high extreme of grid
                    self.H[:,-1,:,0] -= self.chxe[:,-1,:] * ((self.back_extyEz(True) - self.E[:,-1,:,2])/squeeze(self.kappa[...,1]) + self.psiH[:,-1,:,1,0]) #Hx
                    self.H[:,-1,:,2] += self.chze[:,-1,:] * ((self.back_extyEx(True) - self.E[:,-1,:,0])/squeeze(self.kappa[...,1]) + self.psiH[:,-1,:,1,2]) #Hz
            else: #Transversal boundaries
                self.H[:,-1,:,0] -= self.chxe[:,-1,:] * ((self.trans_extyEz(True) - self.E[:,-1,:,2])/squeeze(self.kappa[...,1]) + self.psiH[:,-1,:,1,0]) #Hx
                self.H[:,-1,:,2] += self.chze[:,-1,:] * ((self.trans_extyEx(True) - self.E[:,-1,:,0])/squeeze(self.kappa[...,1]) + self.psiH[:,-1,:,1,2]) #Hz
        
        #Update on z axis
        if self.Nz > 1:
            self.H[:,:,:-1,0] += self.chxe[:,:,:-1] * ((self.E[:,:,1:,1] - self.E[:,:,:-1,1])/squeeze(self.kappa[...,2]) + self.psiH[:,:,:-1,2,0]) #Hx
            self.H[:,:,:-1,1] -= self.chye[:,:,:-1] * ((self.E[:,:,1:,0] - self.E[:,:,:-1,0])/squeeze(self.kappa[...,2]) + self.psiH[:,:,:-1,2,1]) #Hy
            #Boundaries
            if self.orientation == "z": #Normal boundaries
                if self.field == "E": #At low extreme of grid, -1 idx is joined to grid
                    self.H[:,:,-1,0] += self.chxe[:,:,-1] * ((self.grid.E[:,:,0,1] - self.E[:,:,-1,1])/squeeze(self.kappa[...,2]) + self.psiH[:,:,-1,2,0]) #Hx
                    self.H[:,:,-1,1] -= self.chye[:,:,-1] * ((self.grid.E[:,:,0,0] - self.E[:,:,-1,0])/squeeze(self.kappa[...,2]) + self.psiH[:,:,-1,2,1]) #Hy
                elif self.field == "H": #At high extreme of grid
                    self.H[:,:,-1,0] += self.chxe[:,:,-1] * ((self.back_extzEy(True) - self.E[:,:,-1,1])/squeeze(self.kappa[...,2]) + self.psiH[:,:,-1,2,0]) #Hx
                    self.H[:,:,-1,1] -= self.chye[:,:,-1] * ((self.back_extzEx(True) - self.E[:,:,-1,0])/squeeze(self.kappa[...,2]) + self.psiH[:,:,-1,2,1]) #Hy
            else: #Transversal boundaries
                self.H[:,:,-1,0] += self.chxe[:,:,-1] * ((self.trans_extzEy(True) - self.E[:,:,-1,1])/squeeze(self.kappa[...,2]) + self.psiH[:,:,-1,2,0]) #Hx
                self.H[:,:,-1,1] -= self.chye[:,:,-1] * ((self.trans_extzEy(True) - self.E[:,:,-1,0])/squeeze(self.kappa[...,2]) + self.psiH[:,:,-1,2,1]) #Hy
        
        self.previous_H = copy(self.H)
    
    def update_H(self):
        if self.orientation == "x":
            self.grid.H[-1,:,:,1] += self.grid.chye[-1,:,:] * (self.E[0,:,:,2] - self.grid.E[-1,:,:,2]) #Hy
            self.grid.H[-1,:,:,2] -= self.grid.chze[-1,:,:] * (self.E[0,:,:,1] - self.grid.E[-1,:,:,1]) #Hz
        elif self.orientation == "y":
            self.grid.H[:,-1,:,0] -= self.grid.chxe[:,-1,:] * (self.E[:,0,:,2] - self.grid.E[:,-1,:,2]) #Hx
            self.grid.H[:,-1,:,2] += self.grid.chze[:,-1,:] * (self.E[:,0,:,0] - self.grid.E[:,-1,:,0]) #Hz
        elif self.orientation == "z":
            self.grid.H[:,:,-1,0] += self.grid.chxe[:,:,-1] * (self.E[:,:,0,1] - self.grid.E[:,:,-1,1]) #Hx
            self.grid.H[:,:,-1,1] -= self.grid.chye[:,:,-1] * (self.E[:,:,0,0] - self.grid.E[:,:,-1,0]) #Hy
    
    def update_psiE(self):
        self.psiE *= self.b[...,None]
        #Update on x axis
        if self.Nx > 1:
            self.psiE[1:,:,:,0,1] += self.C[1:,:,:,0] * (self.H[1:,:,:,2] - self.H[:-1,:,:,2])
            self.psiE[1:,:,:,0,2] += self.C[1:,:,:,0] * (self.H[1:,:,:,1] - self.H[:-1,:,:,1])
            #Boundaries
            if self.orientation == "x": #Normal boundaries
                if self.field == "E": #At low extreme of grid
                    self.psiE[0,:,:,0,1] += self.C[0,:,:,0] * (self.H[0,:,:,2] - self.back_extxHz(False))
                    self.psiE[0,:,:,0,2] += self.C[0,:,:,0] * (self.H[0,:,:,1] - self.back_extxHy(False))
                elif self.field == "H": #At high extreme of grid, 0 idx is joined to grid
                    self.psiE[0,:,:,0,1] += self.C[0,:,:,0] * (self.H[0,:,:,2] - self.grid.H[-1,:,:,2])
                    self.psiE[0,:,:,0,2] += self.C[0,:,:,0] * (self.H[0,:,:,1] - self.grid.H[-1,:,:,1])
            else: #Transversal boundaries
                self.psiE[0,:,:,0,1] += self.C[0,:,:,0] * (self.H[0,:,:,2] - self.trans_extxHz(False))
                self.psiE[0,:,:,0,2] += self.C[0,:,:,0] * (self.H[0,:,:,1] - self.trans_extxHy(False))
        
        #Update on y axis
        if self.Ny > 1:
            self.psiE[:,1:,:,1,0] += self.C[:,1:,:,1] * (self.H[:,1:,:,2] - self.H[:,:-1,:,2])
            self.psiE[:,1:,:,1,2] += self.C[:,1:,:,1] * (self.H[:,1:,:,0] - self.H[:,:-1,:,0])
            #Boundaries
            if self.orientation == "y": #Normal boundaries
                if self.field == "E": #At low extreme of grid
                    self.psiE[:,0,:,1,0] += self.C[:,0,:,1] * (self.H[:,0,:,2] - self.back_extyHz(False))
                    self.psiE[:,0,:,1,2] += self.C[:,0,:,1] * (self.H[:,0,:,0] - self.back_extyHx(False))
                elif self.field == "H": #At high extreme of grid, 0 idx is joined to grid
                    self.psiE[:,0,:,1,0] += self.C[:,0,:,1] * (self.H[:,0,:,2] - self.grid.H[:,-1,:,2])
                    self.psiE[:,0,:,1,2] += self.C[:,0,:,1] * (self.H[:,0,:,0] - self.grid.H[:,-1,:,0])
            else: #Transversal boundaries
                self.psiE[:,0,:,1,0] += self.C[:,0,:,1] * (self.H[:,0,:,2] - self.trans_extyHz(False))
                self.psiE[:,0,:,1,2] += self.C[:,0,:,1] * (self.H[:,0,:,0] - self.trans_extyHx(False))
        
        #Update on z axis
        if self.Nz > 1:
            self.psiE[:,:,1:,2,0] += self.C[:,:,1:,2] * (self.H[:,:,1:,1] - self.H[:,:,:-1,1])
            self.psiE[:,:,1:,2,1] += self.C[:,:,1:,2] * (self.H[:,:,1:,0] - self.H[:,:,:-1,0])
            #Boundaries
            if self.orientation == "z": #Normal boundaries
                if self.field == "E": #At low extreme of grid
                    self.psiE[:,:,0,2,0] += self.C[:,:,0,2] * (self.H[:,:,0,1] - self.back_extzHy(False))
                    self.psiE[:,:,0,2,1] += self.C[:,:,0,2] * (self.H[:,:,0,0] - self.back_extzHx(False))
                elif self.field == "H": #At high extreme of grid, 0 idx is joined to grid
                    self.psiE[:,:,0,2,0] += self.C[:,:,0,2] * (self.H[:,:,0,1] - self.grid.H[:,:,-1,1])
                    self.psiE[:,:,0,2,1] += self.C[:,:,0,2] * (self.H[:,:,0,0] - self.grid.H[:,:,-1,0])
            else: #Transversal boundaries
                self.psiE[:,:,0,2,0] += self.C[:,:,0,2] * (self.H[:,:,0,1] - self.trans_extzHy(False))
                self.psiE[:,:,0,2,1] += self.C[:,:,0,2] * (self.H[:,:,0,0] - self.trans_extzHx(False))
    
    def update_boundE(self):
        self.update_psiE()
        #Update on x axis
        if self.Nx > 1:
            self.E[1:,:,:,1] -= self.ceyh[1:,:,:] * ((self.H[1:,:,:,2] - self.H[:-1,:,:,2])/squeeze(self.kappa[...,0]) + self.psiE[1:,:,:,0,1]) #Ey
            self.E[1:,:,:,2] += self.cezh[1:,:,:] * ((self.H[1:,:,:,1] - self.H[:-1,:,:,1])/squeeze(self.kappa[...,0]) + self.psiE[1:,:,:,0,2]) #Ez
            #Boundaries
            if self.orientation == "x": #Normal boundaries
                if self.field == "E": #At low extreme of grid
                    self.E[0,:,:,1] -= self.ceyh[0,:,:] * ((self.H[0,:,:,2] - self.back_extxHz(True))/squeeze(self.kappa[...,0]) + self.psiE[0,:,:,0,1]) #Ey
                    self.E[0,:,:,2] += self.cezh[0,:,:] * ((self.H[0,:,:,1] - self.back_extxHy(True))/squeeze(self.kappa[...,0]) + self.psiE[0,:,:,0,2]) #Ez
                elif self.field == "H": #At high extreme of grid, 0 idx is joined to grid
                    self.E[0,:,:,1] -= self.ceyh[0,:,:] * ((self.H[0,:,:,2] - self.grid.H[-1,:,:,2])/squeeze(self.kappa[...,0]) + self.psiE[0,:,:,0,1]) #Ey
                    self.E[0,:,:,2] += self.cezh[0,:,:] * ((self.H[0,:,:,1] - self.grid.H[-1,:,:,1])/squeeze(self.kappa[...,0]) + self.psiE[0,:,:,0,2]) #Ey
            else: #Transversal boundaries
                self.E[0,:,:,1] -= self.ceyh[0,:,:] * ((self.H[0,:,:,2] - self.trans_extxHz(True))/squeeze(self.kappa[...,0]) + self.psiE[0,:,:,0,1]) #Ey
                self.E[0,:,:,2] += self.cezh[0,:,:] * ((self.H[0,:,:,1] - self.trans_extxHy(True))/squeeze(self.kappa[...,0]) + self.psiE[0,:,:,0,2]) #Ez
        
        #Update on y axis
        if self.Ny > 1:
            self.E[:,1:,:,0] += self.cexh[:,1:,:] * ((self.H[:,1:,:,2] - self.H[:,:-1,:,2])/squeeze(self.kappa[...,1]) + self.psiE[:,1:,:,1,0]) #Ex
            self.E[:,1:,:,2] -= self.cezh[:,1:,:] * ((self.H[:,1:,:,0] - self.H[:,:-1,:,0])/squeeze(self.kappa[...,1]) + self.psiE[:,1:,:,1,2]) #Ez
            #Boundaries
            if self.orientation == "y": #Normal boundaries
                if self.field == "E": #At low extreme of grid
                    self.E[:,0,:,0] += self.cexh[:,0,:] * ((self.H[:,0,:,2] - self.back_extyHz(True))/squeeze(self.kappa[...,1]) + self.psiE[:,0,:,1,0]) #Ex
                    self.E[:,0,:,2] -= self.cezh[:,0,:] * ((self.H[:,0,:,0] - self.back_extyHx(True))/squeeze(self.kappa[...,1]) + self.psiE[:,0,:,1,2]) #Ez
                elif self.field == "H": #At high extreme of grid, 0 idx is joined to grid
                    self.E[:,0,:,0] += self.cexh[:,0,:] * ((self.H[:,0,:,2] - self.grid.H[:,-1,:,2])/squeeze(self.kappa[...,1]) + self.psiE[:,0,:,1,0]) #Ex
                    self.E[:,0,:,2] -= self.cezh[:,0,:] * ((self.H[:,0,:,0] - self.grid.H[:,-1,:,0])/squeeze(self.kappa[...,1]) + self.psiE[:,0,:,1,2]) #Ez
            else: #Transversal boundaries
                self.E[:,0,:,0] += self.cexh[:,0,:] * ((self.H[:,0,:,2] - self.trans_extyHz(True))/squeeze(self.kappa[...,1]) + self.psiE[:,0,:,1,0]) #Ex
                self.E[:,0,:,2] -= self.cezh[:,0,:] * ((self.H[:,0,:,0] - self.trans_extyHx(True))/squeeze(self.kappa[...,1]) + self.psiE[:,0,:,1,2]) #Ez
        
        #Update on z axis
        if self.Nz > 1:
            self.E[:,:,1:,0] -= self.cexh[:,:,1:] * ((self.H[:,:,1:,1] - self.H[:,:,:-1,1])/squeeze(self.kappa[...,2]) + self.psiE[:,:,1:,2,0]) #Ex
            self.E[:,:,1:,1] += self.ceyh[:,:,1:] * ((self.H[:,:,1:,0] - self.H[:,:,:-1,0])/squeeze(self.kappa[...,2]) + self.psiE[:,:,1:,2,1]) #Ey
            #Boundaries
            if self.orientation == "z": #Normal boundaries
                if self.field == "E": #At low extreme of grid
                    self.E[:,:,0,0] -= self.cexh[:,:,0] * ((self.H[:,:,0,1] - self.back_extzHy(True))/squeeze(self.kappa[...,2]) + self.psiE[:,:,0,2,0]) #Ex
                    self.E[:,:,0,1] += self.ceyh[:,:,0] * ((self.H[:,:,0,0] - self.back_extzHx(True))/squeeze(self.kappa[...,2]) + self.psiE[:,:,0,2,1]) #Ey
                elif self.field == "H": #At high extreme of grid, 0 idx is joined to grid
                    self.E[:,:,0,0] -= self.cexh[:,:,0] * ((self.H[:,:,0,1] - self.grid.H[:,:,-1,1])/squeeze(self.kappa[...,2]) + self.psiE[:,:,0,2,0]) #Ex
                    self.E[:,:,0,1] += self.ceyh[:,:,0] * ((self.H[:,:,0,0] - self.grid.H[:,:,-1,0])/squeeze(self.kappa[...,2]) + self.psiE[:,:,0,2,1]) #Ey
            else: #Transversal boundaries
                self.E[:,:,0,0] -= self.cexh[:,:,0] * ((self.H[:,:,0,1] - self.trans_extzHy(True))/squeeze(self.kappa[...,2]) + self.psiE[:,:,0,2,0]) #Ex
                self.E[:,:,0,1] += self.ceyh[:,:,0] * ((self.H[:,:,0,0] - self.trans_extzHx(True))/squeeze(self.kappa[...,2]) + self.psiE[:,:,0,2,1]) #Ey
        
        self.previous_E = copy(self.E)
    
    def update_E(self):
        if self.orientation == "x":
            self.grid.E[0,:,:,1] -= self.grid.ceyh[0,:,:] * (self.grid.H[0,:,:,2] - self.H[-1,:,:,2]) #Ey
            self.grid.E[0,:,:,2] += self.grid.cezh[0,:,:] * (self.grid.H[0,:,:,1] - self.H[-1,:,:,1]) #Ez
        elif self.orientation == "y":
            self.grid.E[:,0,:,0] += self.grid.cexh[:,0,:] * (self.grid.H[:,0,:,2] - self.H[:,-1,:,2]) #Ex
            self.grid.E[:,0,:,2] -= self.grid.cezh[:,0,:] * (self.grid.H[:,0,:,0] - self.H[:,-1,:,0]) #Ez
        elif self.orientation == "z":
            self.grid.E[:,:,0,0] -= self.grid.cexh[:,:,0] * (self.grid.H[:,:,0,1] - self.H[:,:,-1,1]) #Ex
            self.grid.E[:,:,0,1] += self.grid.ceyh[:,:,0] * (self.grid.H[:,:,0,0] - self.H[:,:,-1,0]) #Ey
