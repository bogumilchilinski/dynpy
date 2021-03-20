from sympy import *
from sympy.physics.mechanics import *

from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator

from sympy.physics.mechanics import *
from sympy.physics.vector import *

class SDOFHarmonicOscillator(HarmonicOscillator):
    """Ready to use sample Single Degree of Freedom System with mass on spring
    
    """
    def __init__(self, m=Symbol('m',positive=True), k=Symbol('k',positive=True), ivar=Symbol('t'), qs=[z=dynamicsymbols('z')]):
        
        self.m = m
        self.k = k
        self.qs = qs
        
        self.mass = MaterialPoint(m,pos1=qs)
        self.spring = Spring(k,pos1=qs)
        system = self.mass + self.spring
        
        super().__init__(system)