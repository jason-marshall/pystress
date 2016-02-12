import numpy as np
import matplotlib.pyplot as plt
import sys
from enum import Enum
from math import tan, radians

class Dimension(Enum):
    evaluate = 0
    OneD = 1
    TwoD = 2
    ThreeD = 3

class Symmetric(Enum):
    no = 0
    yes = 1
    evaluate = 2

class StressState(object):
    """A stress state in 2D or 3D.

    Attributes:
        stress: Stress state.
        principals: Numpy array of ordered principal stresses.
        eigenvectors: Numpy array of eigenvectors of principal stresses.
        dimension: Dimension of stress state.
        symmetric: Is stress state symmetric?
        max_shear: Maximum shear stress
        hydro_stress: Hydrostatic stress
        dev_stress: Deviatoric stress
        vm_stress: Von-Mises stress criterion
    """

    def __init__(self,
                 s11 = 0.0,
                 s22 = 0.0,
                 s12 = 0.0,
                 s21 = 0.0,
                 s33 = 0.0,
                 s13 = 0.0,
                 s23 = 0.0,
                 s31 = 0.0,
                 s32 = 0.0,
                 dimension = Dimension.evaluate,
                 symmetric = Symmetric.evaluate):
        """Return a valid stress state."""

        # figure out dimension if not given
        if dimension is Dimension.evaluate:
            if s33 == 0.0 and s32 == 0.0 and s31 == 0.0 and s13 == 0.0 and s23 == 0.0:
                self.dimension = Dimension.TwoD
            else:
                self.dimension = Dimension.ThreeD
        elif dimension is Dimension.OneD:
            print "Not much can be done with 1D state, currently not supported!"
            sys.exit()
        else:
            self.dimension = dimension

        # figure out if stress state is symmetric
        if symmetric is Symmetric.evaluate and self.dimension is Dimension.TwoD:
            if s12 == s21:
                self.symmetric = Symmetric.yes
            else:
                self.symmetric = Symmetric.no
        elif symmetric is Symmetric.evaluate and self.dimension is Dimension.ThreeD:
            if s12 == s21 and s13 == s31 and s23 == s32:
                self.symmetric = Symmetric.yes
            else:
                self.symmetric = Symmetric.no
        else:
            self.symmetric = symmetric

        # create stress tensor
        self.stress = np.zeros((self.dimension, self.dimension))
        self.stress[0][0] = s11
        self.stress[0][1] = s12
        self.stress[1][1] = s22
        if self.dimension is Dimension.TwoD:
            if self.symmetric is Symmetric.yes:
                self.stress[1][0] = s12
            else:
                self.stress[1][0] = s21
        else:
            self.stress[0][2] = s13
            self.stress[1][2] = s23
            self.stress[2][2] = s33
            if self.symmetric is Symmetric.yes:
                self.stress[1][0] = s12
                self.stress[2][0] = s13
                self.stress[2][1] = s23
            else:
                self.stress[1][0] = s21
                self.stress[2][0] = s31
                self.stress[2][1] = s32

        # calculate principal stresses and eigenvectors
        evalue,evec = np.linalg.eig(self.stress)
        ev_list = zip( evalue, evec )
        ev_list.sort(key=lambda tup:tup[0], reverse=True)
        self.principals, self.eigenvectors = zip(*ev_list)

        # calculate various properties
        self.hydro_stress = 1.0/3.0*np.sum(self.principals)

        temp_stress = np.zeros((self.dimension, self.dimension))
        np.fill_diagonal(temp_stress,1.0)
        temp_stress = np.multiply(temp_stress, self.hydro_stress)
        self.dev_stress = np.subtract(self.stress, temp_stress)

        if self.dimension is Dimension.TwoD:
            var1 = self.principals[0]-self.principals[1]
            self.vm_stress = np.sqrt(0.5*(var1*var1))
        else:
            var1 = self.principals[0]-self.principals[1]
            var2 = self.principals[1]-self.principals[2]
            var3 = self.principals[2]-self.principals[0]
            self.vm_stress = np.sqrt(0.5*(var1*var1 + var2*var2 + var3*var3))
        
        self.max_shear = (np.max(self.principals) - np.min(self.principals)) / 2.0

    def deviatoric(self):
        """Return the deviatoric stress."""
        return self.dev_stress

    def hydrostatic(self):
        """Return the hydrostatic stress."""
        return self.hydro_stress

    def maxshear(self):
        """Return the max shear stress."""
        return self.max_shear

    def principals(self):
        """Return the principal stresses."""
        return self.principals

    def eigenvectors(self):
        """Return the principal eigenvectors."""
        return self.eigenvectors

    def von_mises(self):
        """Return the von mises stress."""
        return self.vm_stress

    def stress(self):
        """Return the stress."""
        return self.stress

    def dimension(self):
        """Return the dimension."""
        return self.dimension

    def symmetric(self):
        """Return whether the stress tensor is symmetric."""
        return self.symmetric

    def plot_stress(self):
        """Plot the stress state."""
        self.fig = plt.figure(facecolor='w', figsize=(4,3), dpi=600)
        self.ax = plt.subplot(aspect='equal')
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        if self.dimension is Dimension.TwoD:
            circle_center = (np.max(self.principals) + np.min(self.principals)) / 2.0
            circle = plt.Circle((circle_center,0.0),self.max_shear,color='k',fill=False)
            self.fig.gca().add_artist(circle)
        else:
            circle1_center = (self.principals[0] + self.principals[1]) / 2.0
            circle2_center = (self.principals[1] + self.principals[2]) / 2.0
            circle3_center = (self.principals[0] + self.principals[2]) / 2.0
            circle1_radius = (self.principals[0] - self.principals[1]) / 2.0
            circle2_radius = (self.principals[1] - self.principals[2]) / 2.0
            circle3_radius = (self.principals[0] - self.principals[2]) / 2.0
            circle1 = plt.Circle((circle1_center,0.0),circle1_radius,color='k',fill=False)
            circle2 = plt.Circle((circle2_center,0.0),circle2_radius,color='k',fill=False)
            circle3 = plt.Circle((circle3_center,0.0),circle3_radius,color='k',fill=False)
            self.fig.gca().add_artist(circle1)
            self.fig.gca().add_artist(circle2)
            self.fig.gca().add_artist(circle3)

        plt.title("Stress State",fontsize=12, y=1.03)
        plt.xlabel("$\sigma$",fontsize=10)
        plt.ylabel("$\\tau$",fontsize=10)
        self.ax.tick_params(labelsize=10)

        absolute = np.abs(self.principals)
        maxi = np.max(absolute)
        xlim = 0.1*maxi

        if self.principals[-1]-xlim < 0.0:
            xmin = self.principals[-1]-xlim
        else:
            xmin = 0.0
        
        plt.xlim(xmin,self.principals[0]+xlim)
        plt.ylim(0.0,1.1*self.max_shear)
        plt.grid(b=True, which='major', color='gray', ls='--',axis='both',alpha=0.4)
        
        return

    def plot_mohr_coulomb(self, c = 0.0, phi = 30.0):
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        xlims = plt.xlim()
        slope = np.tan(radians(phi))
        p1 = [-c/slope, 0.0]
        p2 = [xlims[1], slope*xlims[1]+c]

        plt.plot([p1[0],p2[0]],[p1[1],p2[1]],'-k',lw=1)



