from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Derivative,Integral, Expr,Function,latex, Heaviside, atan, pi)
from numbers import Number
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex

import numpy as np

from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator, GeometryScene,base_frame,base_origin
from ..utilities.templates import tikz

#from dgeometry import GeometryScene,np

from pylatex import TikZ,TikZNode
import matplotlib.pyplot as plt

import base64
import random
import IPython as IP







class GeometryOfPoint:

    def __init__(self, *args , frame=base_frame, base_origin=base_origin , ivar=Symbol('t')):
        
        #print(type(args),args)
        #name = str(args[0])

        if isinstance(args[0], Point):
            self.Point = args[0]

        elif isinstance(args[0], Expr) or isinstance(args[0], Number):
            P1 = Point('P1')
            P1.set_pos(base_origin, frame.x*args[0] )
            P1.set_vel(frame, frame.x*diff(args[0], ivar) )

            self.Point = P1
            #print(self.Point)

            #print('Expression based on the generalized coordinate')

        elif isinstance(args, tuple):

            if len(args) == 2:
                P1 = Point('P1')
                P1.set_pos(base_origin, frame.x*args[0] + frame.y*args[1])
                P1.set_vel(frame, frame.x*diff(args[0], ivar) + frame.y*diff(args[1], ivar))
                self.Point = P1
                #print(self.Point)

            elif len(args) == 3:
                P2 = Point('P2')
                P2.set_pos(base_origin, frame.x*args[0] + frame.y*args[1] + frame.z*args[2])
                P2.set_vel(frame, frame.x*diff(args[0], ivar) + frame.y*diff(args[1], ivar) + frame.z*diff(args[2], ivar))
                self.Point = P2
                #print(self.Point)

            else:
                print('Tuple out of range')
                P_na = Point('P_na')
                P_na.set_pos(base_origin, frame.x*0)
                P_na.set_vel(frame, frame.x*diff(0, ivar))
                self.Point = P_na

    def get_point(self):
        return self.Point




class Element(LagrangesDynamicSystem):
    """Base class for all elements
    """

    scheme_name = 'damper.PNG'
    real_name = 'damper.PNG'
    _default_folder_path = "./dynpy/models/images/"
    

    
    @classmethod
    def _scheme(cls):

        path = cls._default_folder_path + cls.scheme_name
        
        return path

    @classmethod
    def _real_example(cls):
        path = cls._default_folder_path + cls.real_name

        return path
    
    @classmethod
    def _detail_real(cls):
        path = cls._default_folder_path + cls.detail_real_name

        return path
    
    @classmethod
    def _detail_scheme(cls):
        path = cls._default_folder_path + cls.detail_scheme_name

        return path

    def _init_from_components(self, *args, system=None, **kwargs):

        if system is None:
            composed_system = self._elements_sum
        else:
            composed_system = system

        #print('CS',composed_system._components)
        super().__init__(None, system=composed_system)

        #print('self',self._components)
        if self._components is None:
            comps = {}
        else:
            comps = self._components

        self._components = {**comps, **self.components}

    def _all_default_data(self):

        
        return self.get_default_data()
    
    def _all_numerical_data(self):
        
        return self.get_numerical_data()
        
    
    def get_default_data(self):
        return {}

    def get_numerical_data(self):
        return {}
    
    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict

    def get_numerical_parameters(self):

        default_data_dict = self.get_numerical_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict
    
    @classmethod
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        else:
            path = cls._scheme()

        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))


    def __str__(self):
        

        return f'{(self._label)}'    

    def _tikz_2d(self, language='en',*args,**kwargs):


        
        class_name = self.__class__.__name__


        #print(f'{class_name} element check', self.scheme_options )
        return [TikZNode('Material Point',options=['draw'],text=f'{class_name}')]

    
#Pati
class MaterialPoint(Element):
    """
    Model of a Material point with changing point of mass:

    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'
    
    m=Symbol('M',positive=True)
    pos1=Symbol('pos1',positive=True)
    
    z=dynamicsymbols('z')
    
    def __init__(self,
                 m, 
                 pos1 ,
                 qs=None,
                 frame=base_frame,
                 ivar=Symbol('t'),
                 **kwargs):
        if m is not None: self.m=m
        if pos1 is not None: self.pos1=pos1
        if not qs:
            self.qs = [pos1]

        pos1=GeometryOfPoint(pos1).get_point()
            
        if isinstance(pos1, Point):
            
            v=((base_origin.pos_from(pos1).diff(self.ivar,frame).magnitude())**2)
            
            Lagrangian = S.Half * m * ((base_origin.pos_from(pos1).diff(ivar,frame).magnitude())**2)
            #Lagrangian = S.Half * m * ((base_origin.vel(frame).magnitude())**2)
            
        else:
            Lagrangian = S.Half * m * diff(pos1,ivar)**2


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar,frame=frame,**kwargs)
        
        self._kinetic_energy = Lagrangian
    def symbols_description(self):
        self.sym_desc_dict={
            self.m :r'Mass of body',
            self.pos1 :r'Instantenous position of the body'
        }
        return {**super().symbols_description(),**self.sym_desc_dict}
    
    def get_default_data(self):
        t=self.ivar
        m0=Symbol('m_0',positive=True)
        self.default_data_dict={
            self.m:[m0*no for no in range(1,3)],
        }
        return {**super().get_default_data(),**self.default_data_dict}
    

    def get_numerical_data(self):

        self.default_data_dict={
            self.m:[5* no/100 for no in range(80,120)],
        }
        return self.default_data_dict

    @classmethod
    def from_velocity(cls,
                      m,
                      velocity,
                      qs=None,
                      z=None,
                      ivar=Symbol('t'),
                      frame = base_frame,
                      **kwargs):
        pos1=Integral(velocity,ivar)
        return cls(m=m, pos1=pos1, qs=qs, ivar=ivar, frame = frame)
    @classmethod
    def from_default(cls,
                      m=None,
                      pos1=None,
                      qs=None,
                      ivar=Symbol('t'),
                      frame = base_frame,
                      **kwargs):
        if m is None: m=cls.m
        if pos1 is None: pos1=dynamicsymbols('x')
        if not qs:
            qs=[pos1]
        ivar=cls.ivar
        return cls(m=m, pos1=pos1, qs=qs, ivar=ivar, frame = frame)
    
    def _plot_2d(self, language='en',*args,**kwargs):
      
        class_name =self.__class__.__name__
        
        sch_opt=self.scheme_options
        
        x_shift,y_shift = 0,0
        
        if sch_opt is not None:
            if 'at' in sch_opt:
                x_shift,y_shift = sch_opt['at']
        
        coords = [(0,0),(0,10),(5,10),(5,0),(0,0)]
        x_coords = [  x+x_shift for x,y in  coords ]
        y_coords = [  y+y_shift for x,y in coords]
        
#         my_fun=lambda arg: arg[1]
#         y_coords = list(map(my_fun,  coords))   #[  y for x,y in coords]
#         print(y_coords)
    
        res = GeometryScene.ax_2d.plot(x_coords,y_coords,label=class_name)
        GeometryScene.ax_2d.text(2.5,5,f'${self.m}$',rotation='horizontal',multialignment='center')
        
        

#Mateusz
class Spring(Element):
    """
    Model of a Spring: Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s),
    which analytically display the dynamics of displacing spring after cummulating PE.
    """
    scheme_name = 'spring.PNG'
    real_name = 'spring_real.png'
    
    stiffness=Symbol('k',positive=True)
    l_s=Symbol('l_s', positive=True)
    k0 =Symbol('k_0',positive=True)
    
    def __init__(self, stiffness, pos1, pos2=0, l_s=None ,  qs=None, ivar=Symbol('t'), frame = base_frame,**kwargs):

        if stiffness is not None: self.stiffness=stiffness
        if pos1 is not None: self._pos1=pos1
        self._pos2=pos2
        if l_s is not None: self.l_s=l_s
        else:
            l_s=0
            self.l_s=l_s

        if not qs:
            qs = [pos1]

        Lagrangian = (-S.Half * self.stiffness * ((self.pos2.pos_from(self.pos1)).magnitude() - self.l_s  )**2 )


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, frame=frame,**kwargs)

        self._potential_energy = - Lagrangian
        self._frame = frame

    @property
    def pos1(self):
        return GeometryOfPoint(self._pos1).get_point()

    @property
    def pos2(self):
        return GeometryOfPoint(self._pos2).get_point()

    def force(self):
        if self.pos2==0:
            force = (self.stiffness * self.pos1)
        elif  self.l_s == 0:
            force = self.stiffness * (self.pos2.pos_from(self.pos1)) & self._frame.x
            
        else:
            force = self.stiffness * (self.pos2.pos_from(self.pos1)).magnitude() - self.l_s
        return force
    

    def _plot_2d(self, language='en',*args,**kwargs):

        class_name = self.__class__.__name__

        coords=[(0,0),(0,3),(3,6),(-3,12),(0,15),(0,18)]

        sch_opt=self.scheme_options
        x_shift,y_shift = 0,0
        
        if sch_opt is not None:
            if 'at' in sch_opt:
                x_shift,y_shift = sch_opt['at']

        x_span = np.arange(0,21,1)


        x_coords = np.sin(x_span*np.pi/2)+x_shift
        y_coords = (x_span*0.8)+y_shift
        
        x_end=x_coords[-1]
        y_end=y_coords[-1]

        #print(y_coords)
        
        res = GeometryScene.ax_2d.plot(x_coords,y_coords,label=class_name)
        res =GeometryScene.ax_2d.text(x_end,y_end,'$' + latex(self.stiffness) + '$')

        
        
    def get_default_data(self):

        k0= self.k0
        self.default_data_dict={
            self.stiffness: [S.One *  k0 * no for no in range(80, 135)],
        }
        return {**super().get_default_data(),**self.default_data_dict}
        
    def get_numerical_data(self):

        self.default_numerical_dict={
            self.stiffness: [S.One / 100 *  50e3 * no for no in range(80, 135)],
        }
        
        return {**super().get_numerical_data(),**self.default_numerical_dict}
    
    @classmethod 
    def from_default(cls,
                 stiffness=None,
                 pos1=None,
                 frame=base_frame,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if stiffness is not None: cls.stiffness=stiffness
        if pos1 is None: pos1= dynamicsymbols('z') 
        if not qs: 
            qs = [pos1] 
        ivar=cls.ivar

        return cls(stiffness=stiffness, pos1=pos1, qs=qs,ivar=ivar, frame = frame)

class EngineMount(Spring):

    """
    Model of a Engine Mount: Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s),
    which analytically display the dynamics of displacing spring after cummulating PE.
    """
    real_path='buick_regal_3800.jpg'#buick_regal_3800.jpg
    scheme_name = 'spring.PNG'
    def pin_diameter(self):
        #średnica sworznia ze wzoru na ścinanie
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)

        return sqrt((4*S.One*self.force())/(pi*kt*Re)).doit()

    def pin_force(self):

        return self.force()

    def get_numerical_data(self):
        self.default_data_dict={
            self.stiffness: [1e5 * no for no in range(1,8)]
        }
        return self.default_data_dict
    
    def _plot_2d(self, language='en',*args,**kwargs):

        class_name = self.__class__.__name__

        coords=[(0,0),(0,3),(3,6),(-3,12),(0,15),(0,18)]

        sch_opt=self.scheme_options
        x_shift,y_shift = 0,0
        
        if sch_opt is not None:
            if 'at' in sch_opt:
                x_shift,y_shift = sch_opt['at']

        x_span = np.arange(0,21,1)


        x_coords = np.sin(x_span*np.pi/2)+x_shift
        y_coords = (x_span*0.4)+y_shift
        
        x_end=x_coords[0]+2
        y_end=y_coords[0]/2

        #print(y_coords)
        
        res = GeometryScene.ax_2d.plot(x_coords,y_coords,label=class_name)
        res = GeometryScene.ax_2d.plot(x_end,y_end,'*r')
        res =GeometryScene.ax_2d.text(x_end,y_end,"pin")
        res =GeometryScene.ax_2d.text(x_end,y_end,'$' + latex(self.stiffness) + '$')
    
class DampedEngineMount(EngineMount):

    """
    Model of a Damped Engine Mount: Creates a singular model, after inputing correct values of stiffeness - k, damping coefficient - c and general coordinate(s),
    which analytically display the dynamics of displacing spring after cummulating PE and damper after commulating KE.
    """
    
    real_path='buick_regal_3800.jpg' #buick_regal_3800.jpg
    scheme_name = 'spring.PNG'
    scheme_name = 'damped_engine_vertical_spring_fy.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'
    
    stiffness = Symbol('k', positive=True)
    damping_coeff = Symbol('c', positive=True)
    _pos1 = dynamicsymbols('z')
    #_pos2 = Symbol('pos_2', positive=True)
    _pos2 = 0
    l_s = Symbol('l_s', positive= True)
    ivar = Symbol('t')
    
    def __init__(
            self,
            stiffness,
            damping_coeff,
            pos1,
            pos2=0,
            l_s=None,
            qs=None,
            ivar=Symbol('t'),
            **kwargs):

        if stiffness is not None: self.stiffness = stiffness
        if damping_coeff is not None:  self.damping_coeff = damping_coeff
        if pos1 is not None: self._pos1=pos1
        if pos2 == 0: self._pos2 = 0
        else: self._pos2 = pos2
        if qs is not None: self.qs=qs
        else: self.qs=[self._pos1]
        if l_s is not None: self.l_s=l_s
        else:
            l_s=0
            self.l_s=l_s
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        
        self._engine_mount = EngineMount(self.stiffness, pos1=self.pos1, pos2=self.pos2, l_s=self.l_s, qs=self.qs)(label='Elasticity',scheme_options={'at':(2,4)})
        self._damper = Damper(self.damping_coeff, pos1=self.pos1, pos2=self.pos2, qs=self.qs)(label='Damping',scheme_options={'at':(-2,4)})

        components['_engine_mount'] = self._engine_mount
        components['_damper'] = self._damper

        return components

    def force(self):
        return self._engine_mount.force() + self._damper.force()         #WorkInProgress

    @classmethod 
    def from_default(cls,
                stiffness=None,
                damping_coeff=None,
                pos1=None,
                pos2=0,
                l_s=None,
                qs=None,
                ivar=Symbol('t'),
                 **kwargs):

        return cls(stiffness=stiffness, damping_coeff=damping_coeff, pos1=pos1, pos2=pos2, l_s=l_s, qs=qs,ivar=ivar)

    def _plot_2d(self, language='en',*args,**kwargs):

        class_name = self.__class__.__name__

        coords=[(0,0),(0,3),(3,6),(-3,12),(0,15),(0,18)]

        sch_opt=self.scheme_options

        x_shift,y_shift = 0,0

        if sch_opt is not None:
            if 'at' in sch_opt:
                x_shift,y_shift = sch_opt['at']

        x_span = np.arange(0,21,1)

        x_coords = np.sin(x_span*np.pi/2)+x_shift
        y_coords = (x_span*0.4)+y_shift
        
        x_end=x_coords[0]+2
        y_end=y_coords[0]/2

        #print(y_coords)
#         for comp in self.components.values():
#             comp._plot_2d(language=language, *args, **kwargs)
        res = GeometryScene.ax_2d.plot(x_coords,y_coords,label=class_name)
        res = GeometryScene.ax_2d.plot(x_end,y_end,'*r')
        res =GeometryScene.ax_2d.text(x_end,y_end,"pin")
        res =GeometryScene.ax_2d.text(x_end,y_end,'$' + latex(self.stiffness) + '$')
    
    
class GravitationalForce(Element):
    """
    Model of a changing centroid for potential energy. Creates a singular model, after inputing correct values of gravity field - g,
     mass of - m as well as additionaly the general coordiante
    """
    scheme_name = ''
    real_name = ''
    def __init__(self, m, g, pos1, qs=None, ivar=Symbol('t')):

        if not qs:
            qs=[pos1]
        self.g = g

        Lagrangian = -(m * g * (pos1))
        


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        self._potential_energy = - Lagrangian
        
    def get_numerical_data(self):
        self.default_data_dict={
            self.g: [9.81]
        }
        return self.default_data_dict


class Disk(Element):
    """
    Model of a Disk:
    Creates a singular model, after inputing correct values of moment of inertia - I and rotational general coordinate, which analytically displays the dynamics of a rotating wheel.
    """
    scheme_name = 'disk.png'
    real_name = 'disk.png'
    def __init__(self, I, pos1 , qs=None,  frame=base_frame, ivar=Symbol('t')):
        
        if not qs:
            self.qs = [pos1]

        pos1=GeometryOfPoint(pos1).get_point()

        if isinstance(pos1, Point):
            Lagrangian = S.Half * I * ((base_origin.pos_from(pos1).diff(ivar,frame).magnitude())**2)
                
        else:
            Lagrangian = S.Half * I * diff(pos1,ivar)**2


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, frame=frame)
        self._kinetic_energy = Lagrangian

        

class RigidBody2D(Element):
    """
    Model of a 2DoF Rigid body:
    """
    """
    Creates a singular model, after inputing correct values of moment of inertia - I for rotational element and mass - m for linear element and general coordinates, which             analytically display the dynamics of rotating and propeling object: beam.
    """
    scheme_name = 'rigid_body2D.png'
    real_name = 'rigid_body2D.png'
    def __init__(self, m, I, pos_lin=0, pos_rot=0, qs=None, frame=base_frame, ivar=Symbol('t')):
        
        if pos_lin == 0 and pos_rot == 0:
            self.qs = [pos_lin, pos_rot]

        if isinstance(pos_lin, Point) and isinstance(pos_rot, Point):
            Lagrangian = S.Half * m * ((base_origin.pos_from(pos_lin).diff(ivar,frame).magnitude())**2) + S.Half * I *((base_origin.pos_from(pos_rot).diff(ivar,frame).magnitude())**2)
        
        else:
            Lagrangian = S.Half * m * diff(pos_lin,ivar)**2 + S.Half * I * diff(pos_rot,ivar)**2

        
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, frame=frame)
        
        self._kinetic_energy = Lagrangian
        
class Damper(Element):
    """
    Model of a Damper:

    Creates a singular model, after inputing correct values of the damping coefficient - c and general coordinates, which establishes a damping force directly proportional to the velocity, in time, of an inertial virbating element.
    """
    scheme_name = 'damper.PNG'
    real_name = 'damper.PNG'
    
    c= Symbol('c', positive=True)
    
    def __init__(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=base_frame):
        
        if c is not None: self.c=c
        self.pos1=pos1
        self.ivar=ivar
        
        if not qs:
            self.qs = [pos1]

        if isinstance(pos1, Point):
            pos1=GeometryOfPoint(pos1,frame=frame).get_point()
            pos2=GeometryOfPoint(pos2,frame=frame).get_point()
            
            
            
            D = S.Half * c * ((pos1.vel(frame)-pos2.vel(frame)).magnitude()**2).doit()  # pos1.vel(frame).magntidue() behaves as diff(pos1,ivar)
        else:
            D = S.Half * c * (diff(pos1-pos2,ivar))**2

        points_dict={}
        for coord in qs:
            
            coord_vel=diff(coord,ivar)
            
            P_tmp = Point(f'P_{str(coord)}')
            P_tmp.set_vel(frame, coord_vel * frame.x)
            points_dict[coord_vel]=P_tmp

    
        forcelist = [ (point_dmp,-diff(D,coord_vel)*frame.x)  for coord_vel,point_dmp in points_dict.items() ]
        #print(forcelist)
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        self._dissipative_potential = D
        
    def force(self):
        
        return self.c * diff(self.q[0],self.ivar)
     
    def _plot_2d(self, language='en',*args,**kwargs):

        class_name = self.__class__.__name__

        sch_opt=self.scheme_options
        
        comps = self.components
        
        x_shift,y_shift = 0,0
        
        if sch_opt is not None:
            if 'at' in sch_opt:
                x_shift,y_shift = sch_opt['at']

        res = GeometryScene.ax_2d.plot(np.array([0,0,1,-1]) + x_shift+2,np.array([-1,2,2,2]) +  y_shift, label=class_name, color = 'k') + GeometryScene.ax_2d.plot(np.array([1.5,1.5,0,0,0,-1.5,-1.5]) + x_shift+2,np.array([1,3,3,5,3,3,1]) + y_shift, label=class_name, color = 'k')
        
        text = GeometryScene.ax_2d.text(np.array([1.5]),np.array([3]),f'{class_name}',multialignment='center')
        
class Excitation(Element):
    """
    Model of a harmonic extorsion applied onto the elemnt:
    """
    scheme_name = ''
    real_name = ''
    def __init__(self, f, pos_rot, ivar=Symbol('t'), frame=base_frame):
        
        qs = [pos_rot]
        
        pos_rot = GeometryOfPoint(pos_rot).get_point()
        
        dpos_rot = diff(pos_rot,ivar)
        
        P = Point('P')
        P.set_vel(frame, dpos_rot * frame.x)
        P.vel(frame)

        F = f * cos(dpos_rot*ivar) + f * sin(dpos_rot*ivar)
        
        forcelist = [(P, F*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

class Force(Element):
    """
    Creates enforcement.
    """
    scheme_name = ''
    real_name = ''
    def __init__(self,
                 force,
                 pos1 = None,
                 qs = None,
                 ivar=Symbol('t'),
                 frame = base_frame):

        if qs is None:
            qs = [pos1]
            
        else:
            qs = qs
    
        if isinstance(pos1, Point):
            P = pos1
            if qs is None:
                diffs=P.vel(frame).magnitude().atoms(Derivative)
                
                qs= [deriv.args[0] for deriv in diffs]
                print(qs)
        else:
            
            P = Point('P')
            P.set_vel(frame,pos1.diff(ivar)*frame.x)
            force=+force*frame.x

        forcelist=[(P,force)]
        
        #display(forcelist)
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

# class CombustionEngine(Force):

#     eta= Symbol('eta', positive = True)
#     n_engine = Symbol('n')
#     omega=Symbol('omega')

#     def __init__(self, omega, characteristic=None, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):

#         if characteristic is None:
#             self._characteristic= self.characteristic
#         else:
#             self._characteristic =  characteristic
#         if eta is not None: self.eta = eta
#         if omega is not None: self.omega = omega
            

#         super().__init__(self.eta*self._characteristic, omega.integrate(ivar), qs=qs, ivar=ivar)

#     @property
#     def characteristic(self):
#         n_engine = self.omega*60/(2*np.pi)
#         return -4e-5*self.n_engine**2 + 0.2157*self.n_engine + 113.07

#     @classmethod
#     def from_data(cls, data, n, degree=5, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):

#         omega = 2*3.14/60*n
        
#         x = data.index.to_numpy()
#         y = data.iloc[:,0].to_numpy()

#         coeffs=reversed((np.polyfit(x, y, degree)))

        
#         omega_min = 2*3.14/60*x[0]
#         omega_max = 2*3.14/60*x[-1]
        
#         char_slope = 10
        
#         window   =   ((S.One/2+atan(char_slope*(omega-omega_min))/np.pi))    -  ((S.One/2+atan(char_slope*(omega-omega_max))/np.pi))
#         #window = Heaviside(omega-omega_min) - Heaviside(omega-omega_max)
        
#         omega = 2*3.14/60*n
        
#         char_poly  =window *sum([ coeff * n ** no for no, coeff in enumerate(coeffs)])
        

        
#         obj = cls(omega = omega, characteristic = char_poly, eta=eta, qs = qs, ivar = ivar, frame = frame)
#         obj._char_poly  =  char_poly
        
#         return obj

#     @classmethod
#     def from_rpm(cls, n, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):
#         omega = 2*3.14/60*n
#         return cls(omega = omega, eta = eta, qs=qs, ivar=ivar, frame = frame)

class CombustionEngine(Force):

    eta= Symbol('eta', positive = True)
    n_engine = Symbol('n')
    omega=Symbol('omega')
    n_min = Symbol('n_min')
    n_max = Symbol('n_max')

    def __init__(self, omega, characteristic=None, eta = None, n_min=None , n_max=None, qs=None, ivar=Symbol('t'),frame = base_frame):

        if characteristic is None:
            self._characteristic= self.characteristic
        else:
            self._characteristic =  characteristic
        if eta is not None: self.eta = eta
        if omega is not None: self.omega = omega
        if n_min is not None: self.n_min = n_min
        if n_max is not None: self.n_max = n_max

        super().__init__(self.eta*self._characteristic, omega.integrate(ivar), qs=qs, ivar=ivar)

    @property
    def characteristic(self):
        n_engine = self.omega*60/(2*np.pi)
        return -4e-5*self.n_engine**2 + 0.2157*self.n_engine + 113.07

    @classmethod
    def from_data(cls, data, n, degree=5, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):

        omega = 2*3.14/60*n
        
        x = data.index.to_numpy()
        y = data.iloc[:,0].to_numpy()

        coeffs=reversed((np.polyfit(x, y, degree)))

        
        omega_min = 2*3.14/60*x[0]
        omega_max = 2*3.14/60*x[-1]
        
        char_slope = 10e60
        
        window   =   ((S.One/2+atan(char_slope*(omega-omega_min))/np.pi))    -  ((S.One/2+atan(char_slope*(omega-omega_max))/np.pi))
        #window = Heaviside(omega-omega_min) - Heaviside(omega-omega_max)
        
        omega = 2*3.14/60*n
        
        char_poly  =window *sum([ coeff * n ** no for no, coeff in enumerate(coeffs)])
        

        
        obj = cls(omega = omega, characteristic = char_poly, eta=eta, qs = qs, ivar = ivar, frame = frame)
        obj._char_poly  =  char_poly
        
        return obj
    
    @classmethod
    def from_data_gearbox(cls, data, n, degree=5, eta = None, n_min=None , n_max=None , qs=None, ivar=Symbol('t'),frame = base_frame):

        omega = 2*3.14/60*n
        
        x = data.index.to_numpy()
        y = data.iloc[:,0].to_numpy()

        coeffs=reversed((np.polyfit(x, y, degree)))

        if n_min is not None:
            omega_min =  2*3.14/60*n_min
        else:
            omega_min = 2*3.14/60*x[0]
        if n_min is not None:
            omega_max =  2*3.14/60*n_max
        else:
            omega_max = 2*3.14/60*x[0]
        
        char_slope = 10e60
        
        window   =   ((S.One/2+atan(char_slope*(omega-omega_min))/np.pi))    -  ((S.One/2+atan(char_slope*(omega-omega_max))/np.pi))
        #window = Heaviside(omega-omega_min) - Heaviside(omega-omega_max)
        
        omega = 2*3.14/60*n
        
        char_poly  =window *sum([ coeff * n ** no for no, coeff in enumerate(coeffs)])
        
        obj = cls(omega = omega, characteristic = char_poly, eta=eta, qs = qs, ivar = ivar, frame = frame)
        obj._char_poly  =  char_poly
        
        return obj
    
    @classmethod
    def from_data_raw(cls, data, n, degree=5, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):

        omega = 2*3.14/60*n
        
        x = data.index.to_numpy()
        y = data.iloc[:,0].to_numpy()

        coeffs=reversed((np.polyfit(x, y, degree)))

        
        omega_min = 2*3.14/60*x[0]
        omega_max = 2*3.14/60*x[-1]
        
        char_slope = 10e60
        
        window = 1
        #window   =   ((S.One/2+atan(char_slope*(omega-omega_min))/np.pi))    -  ((S.One/2+atan(char_slope*(omega-omega_max))/np.pi))
        #window = Heaviside(omega-omega_min) - Heaviside(omega-omega_max)
        
        omega = 2*3.14/60*n
        
        char_poly  =window *sum([ coeff * n ** no for no, coeff in enumerate(coeffs)])
        

        
        obj = cls(omega = omega, characteristic = char_poly, eta=eta, qs = qs, ivar = ivar, frame = frame)
        obj._char_poly  =  char_poly
        
        return obj
    
    @classmethod
    def from_rpm(cls, n, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):
        omega = 2*3.14/60*n
        return cls(omega = omega, eta = eta, qs=qs, ivar=ivar, frame = frame)
    
class IntegralElement(Force):
    """
    Model of an Integral Element :
        Creates a additional force which allows to control output value.

     Examples
        =======
        An implementation of the control system's integral element:
        
        Example No.1 - error is defined as the expresion:
        >>> t = symbols('t')
        >>> I,reference = symbols('I reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> error = reference - x
        >>> IntegralElement(I , error = error , target = x , qs=[x])
        
        Example No.2 - direct parameter regulation:
        >>> t = symbols('t')
        >>> I, reference= symbols('I, reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> IntegralElement(I , error = x  , reference = reference , qs=[x])
        
        Example No.3 - indirect parameter regulation:
        >>> t = symbols('t')
        >>> I, reference= symbols('I, reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> vel = diff(x,t)
        >>> IntegralElement(I , error = vel , target = x , reference = reference, qs=[x])

        -We define the symbols and dynamicsymbols
        -We choose regulation strategy
        -We initiate an object IntegralElement in proper way
    
    """
    k_I = Symbol('k_I', positive=True)
    error = dynamicsymbols('e')
    target = dynamicsymbols('target')
    reference = S.Zero

    def __init__(self, k_I, error = None, target = None, reference = None,  qs=None, ivar=Symbol('t'), frame = base_frame):

        if reference is not None: self.reference = reference
        if error is not None: self.error = error

        if target is None and isinstance(self.error,Function):
               target=self.error

        if isinstance(self.error,Expr) and self.reference == S.Zero:
            super().__init__(k_I*self.error.integrate(ivar), target, qs=qs, ivar=ivar)
        else:
            super().__init__(k_I*(self.reference - self.error).integrate(ivar), target, qs=qs, ivar=ivar)


    def _plot_2d(self, language='en'):

        class_name = self.__class__.__name__

        x_coords = np.array([0,0,5,5,0])
        y_coords = np.array([0,3,3,0,0])

        res = GeometryScene.ax_2d.plot(x_coords,y_coords, label=class_name, color = 'k') + GeometryScene.ax_2d.plot([-1,0],[1.5,1.5], label=class_name, color = 'k') + GeometryScene.ax_2d.plot([6,5],[1.5,1.5], label=class_name, color = 'k')
        
        text = GeometryScene.ax_2d.text(np.array([1]),np.array([1.5]),f'{class_name}',multialignment='center')
        
    def _tikz_scheme(self, language='en'):

        class_name = self.__class__.__name__

        schm = tikz.TikzStandalone()
        tk_node = TikZ()
        
        with schm.create(TikZ()) as tkz:
            with tkz.create(TikZNode('Integral Element', options=['draw'] , text = f'{class_name}')) as node:
                pass
            with tkz.create(TikZNode('Integral Element', options=['draw'] , text = f'{class_name}')) as node:
                pass
        return schm
        
class DerivativeElement(Force):
    """
    Model of an Derivative Element :
        Creates a additional force which allows to control output value.

     Examples
        =======
        An implementation of the control system's derivative element:
        
        Example No.1 - error is defined as the expresion:
        >>> t = symbols('t')
        >>> I,reference = symbols('I reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> error = reference - x
        >>> DerivativeElement(I , error = error , target = x , qs=[x])
        
        Example No.2 - direct parameter regulation:
        >>> t = symbols('t')
        >>> I, reference= symbols('I, reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> DerivativeElement(I , error = x  , reference = reference , qs=[x])
        
        Example No.3 - indirect parameter regulation:
        >>> t = symbols('t')
        >>> I, reference= symbols('I, reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> vel = diff(x,t)
        >>> DerivativeElement(I , error = vel , target = x , reference = reference, qs=[x])

        -We define the symbols and dynamicsymbols
        -We choose regulation strategy
        -We initiate an object DerivativeElement in proper way
    
    """
    k_D = Symbol('k_D', positive=True)
    error = dynamicsymbols('e')
    target = dynamicsymbols('target')
    reference = S.Zero
    
    def __init__(self, k_D, error = None, target=None, reference = None,  qs=None, ivar=Symbol('t'), frame = base_frame):
        
        if reference is not None: self.reference = reference
        if error is not None: self.error = error

        if target is None and isinstance(self.error,Function):
               target=self.error

        if isinstance(self.error,Expr) and self.reference == S.Zero:
            super().__init__(k_D*diff(self.error,ivar), target, qs=qs, ivar=ivar)
        else:
            super().__init__(k_D*diff((self.reference - self.error),ivar), target, qs=qs, ivar=ivar)

    def _plot_2d(self, language='en'):

        class_name = self.__class__.__name__

        x_coords = np.array([0,0,5,5,0])
        y_coords = np.array([0,3,3,0,0])

        res = GeometryScene.ax_2d.plot(x_coords,y_coords, label=class_name, color = 'k') + GeometryScene.ax_2d.plot([-1,0],[1.5,1.5], label=class_name, color = 'k') + GeometryScene.ax_2d.plot([6,5],[1.5,1.5], label=class_name, color = 'k')
        
        text = GeometryScene.ax_2d.text(np.array([1.2]),np.array([1.5]),f'{class_name}',multialignment='center')
        
    def _tikz_scheme(self, language='en'):

        class_name = self.__class__.__name__

        schm = tikz.TikzStandalone()
        tk_node = TikZ()
        
        with schm.create(TikZ()) as tkz:
            with tkz.create(TikZNode('Derivative Element', options=['draw'] , text = f'{class_name}')) as node:
                pass
            with tkz.create(TikZNode('Derivative Element', options=['draw'] , text = f'{class_name}')) as node:
                pass
        return schm


class ProportionalElement(Force):
    """
    Model of an Proportional Element :
        Creates a additional force which allows to control output value.

     Examples
        =======
        An implementation of the control system's proportional element:
        
        Example No.1 - error is defined as the expresion:
        >>> t = symbols('t')
        >>> I,reference = symbols('I reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> error = reference - x
        >>> ProportionalElement(I , error = error , target = x , qs=[x])
        
        Example No.2 - direct parameter regulation:
        >>> t = symbols('t')
        >>> I, reference= symbols('I, reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> ProportionalElement(I , error = x  , reference = reference , qs=[x])
        
        Example No.3 - indirect parameter regulation:
        >>> t = symbols('t')
        >>> I, reference= symbols('I, reference')
        >>> x = dynamicsymbols('x') # Generalized Coordinates
        >>> vel = diff(x,t)
        >>> ProportionalElement(I , error = vel , target = x , reference = reference, qs=[x])

        -We define the symbols and dynamicsymbols
        -We choose regulation strategy
        -We initiate an object ProportionalElement in proper way
    
    """
    k_P = Symbol('k_P', positive=True)
    error = Symbol('e', positive=True)
    target = Symbol('target', positive=True)
    reference = S.Zero
    
    def __init__(self, k_P, error = None, target=None, reference = None,  qs=None, ivar=Symbol('t'), frame = base_frame):

        if reference is not None: self.reference = reference
        if error is not None: self.error = error

        if target is None and isinstance(self.error,Function):
               target=self.error

        if isinstance(self.error,Expr) and self.reference == S.Zero:
            super().__init__(k_P*(self.error), target, qs=qs, ivar=ivar)
        else:
            super().__init__(k_P*(self.reference - self.error), target, qs=qs, ivar=ivar)

    def _plot_2d(self, language='en'):

        class_name = self.__class__.__name__

        x_coords = np.array([0,0,5,5,0])
        y_coords = np.array([0,3,3,0,0])

        res = GeometryScene.ax_2d.plot(x_coords,y_coords, label=class_name, color = 'k') + GeometryScene.ax_2d.plot([-1,0],[1.5,1.5], label=class_name, color = 'k') + GeometryScene.ax_2d.plot([6,5],[1.5,1.5], label=class_name, color = 'k')
        
        text = GeometryScene.ax_2d.text(np.array([1]),np.array([1.5]),f'{class_name}',multialignment='center')
        
    def _tikz_scheme(self, language='en'):

        class_name = self.__class__.__name__

        schm = tikz.TikzStandalone()
        tk_node = TikZ()
        
        with schm.create(TikZ()) as tkz:
            with tkz.create(TikZNode('Proportional Element', options=['draw'] , text = f'{class_name}')) as node:
                pass
            with tkz.create(TikZNode('Proportional Element', options=['draw'] , text = f'{class_name}')) as node:
                pass
        return schm

        
        
class PID(Element):
    """
    Model of a PID controller:

Creates a model of a PID controller (proportional , integral , derivative) which initates a smoother control over delivered system oscillatory system in form of Lagrange method. Taking into account that a PID regulator is ment for a SISO system the model is built for sDoF system, hence in case of building a mDoF the process will require using two PIDs.
    """
    scheme_name = ''
    real_name = ''
    def __init__(self, kp, ki, kd, pos1, qs=None, ivar=Symbol('t'), frame=base_frame):

        if qs:
            qs = [pos1]

        dpos1 = diff(pos1, ivar)

        In = pos1  # In - error
        u = pos1*ki + dpos1*kp + d2pos1*kd
        G = In - u
        
        forcelist = [(P, G*frame.x)]
        
        pos1=GeometryOfPoint(pos1,frame=frame).get_point()
        
        if isinstance(pos1, Point):
            D = S.Half * c * ((pos1.vel(frame)-pos2.vel(frame)).magnitude()**2).doit()  #pos1*ki + pos1.vel(frame)*kp + d2pos1*kd
        else:
            D = S.Half * c * (diff(pos1-pos2,ivar))**2

        points_dict={}
        for coord in qs:
            
            coord_vel=diff(coord,ivar)
            
            P_tmp = Point(f'P_{str(coord)}')
            P_tmp.set_vel(frame, coord_vel * frame.x)
            points_dict[coord_vel]=P_tmp

        forcelist = [ (point_dmp,-diff(G,coord_vel)*frame.x)  for coord_vel,point_dmp in points_dict.items() ]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)