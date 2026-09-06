from functools import cached_property

from sympy import Function, Rational, Symbol, cos, diff, sin
from sympy.physics import units
from sympy.physics.mechanics import dynamicsymbols

from ..elements import Element, Force
from .principles import ComposedSystem

ureg = units


class MasslessElasticShaft(Element):

	scheme_name = "MasslessElasticShaft.png"  # ToBeAdded
	real_name = "MasslessElasticShaft.png"  # ToBeAdded

	def __init__(
		self,
		stiffness=None,
		horizontal_disp=None,
		vertical_disp=None,
		qs=None,
		ivar=Symbol("t"),
		**kwargs
	):
		self.k = stiffness if stiffness is not None else Symbol("k", positive=True)
		self.h = (
			horizontal_disp if horizontal_disp is not None else dynamicsymbols("h")
		)
		self.v = vertical_disp if vertical_disp is not None else dynamicsymbols("v")
		self.ivar = ivar

		if qs is None:
			qs = [self.h, self.v]

		Lagrangian = -Rational(1, 2) * self.k * (
			self.h**2 + self.v**2
		)

		super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, **kwargs)

		self._potential_energy = -Lagrangian


class Rotor(Element):

	def __init__(self, m, I, e, h, v, phi, qs=None, ivar=Symbol("t"), **kwargs):
		self.ivar = ivar

		# Position of the center of mass C
		xc = h + e * cos(phi)
		yc = v + e * sin(phi)

		# Velocity of the center of mass C
		v_xc = diff(xc, ivar)
		v_yc = diff(yc, ivar)

		v_squared = v_xc**2 + v_yc**2

		# Linear kinetic energy
		T_trans = Rational(1, 2) * m * v_squared

		# Rotational kinetic energy
		phi_dot = diff(phi, ivar)
		T_rot = Rational(1, 2) * I * phi_dot**2

		# Total kinetic energy
		Lagrangian = T_trans + T_rot

		super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, **kwargs)


class FreeEngineMDOF(ComposedSystem):
	"""
	Three degrees of freedom (3DOF) rotor model: h, v, phi.
	"""

	scheme_name = "engine_block_XXX.png"  # ToBeAdded
	real_name = "engine_real_XXX.PNG"  # ToBeAdded

	def __init__(
		self,
		m=None,
		I_engine=None,
		e=None,
		k=None,
		h=None,
		v=None,
		phi=None,
		ivar=Symbol("t"),
		**kwargs
	):

		self.m = m if m is not None else Symbol("m", positive=True)
		self.I_engine = I_engine if I_engine is not None else Symbol("I", positive=True)
		self.e = e if e is not None else Symbol("e", positive=True)
		self.k = k if k is not None else Symbol("k", positive=True)

		# Generalized coordinates
		self.h = h if h is not None else dynamicsymbols("h")
		self.v = v if v is not None else dynamicsymbols("v")
		self.phi = phi if phi is not None else dynamicsymbols("\\varphi")

		self.ivar = ivar

		# Generalized coordinates
		self.qs = [self.h, self.v, self.phi]

		self._init_from_components(**kwargs)

	@property
	def components(self):
		components = {}

		self._rotor = Rotor(
			m=self.m,
			I=self.I_engine,
			e=self.e,
			h=self.h,
			v=self.v,
			phi=self.phi,
			qs=self.qs,
			ivar=self.ivar,
		)

		self._shaft = MasslessElasticShaft(
			stiffness=self.k,
			horizontal_disp=self.h,
			vertical_disp=self.v,
			qs=self.qs,
			ivar=self.ivar,
		)

		components["_rotor"] = self._rotor
		components["_shaft"] = self._shaft

		return components

	def symbols_description(self):
		return {
			self.m: r"Mass of the rotor",
			self.I_engine: r"Mass moment of inertia",
			self.e: r"Eccentricity",
			self.k: r"Shaft stiffness",
			self.h: r"Horizontal displacement",
			self.v: r"Vertical displacement",
			self.phi: r"Rotation angle",
		}


class VibratingRotor(ComposedSystem):
	"""
	Rotor model with external torque M(t).
	"""

	scheme_name = "VibratingRotorScheme.png"

	def __init__(
		self,
		m=None,
		I_engine=None,
		e=None,
		k=None,
		M_t=None,
		h=None,
		v=None,
		phi=None,
		qs=None,
		ivar=Symbol("t"),
		**kwargs
	):
		self.ivar = ivar

		if M_t is None:
			self.M_t = Function("M_{ext}")(self.ivar)
		else:
			self.M_t = M_t

		self.m = m if m is not None else Symbol("m", positive=True)
		self.I_engine = I_engine if I_engine is not None else Symbol("I", positive=True)
		self.e = e if e is not None else Symbol("e", positive=True)
		self.k = k if k is not None else Symbol("k", positive=True)

		# Generalized coordinates
		self.h = h if h is not None else dynamicsymbols("h")
		self.v = v if v is not None else dynamicsymbols("v")
		self.phi = phi if phi is not None else dynamicsymbols("\\varphi")

		self.qs = qs if qs is not None else [self.h, self.v, self.phi]

		self._init_from_components(**kwargs)

	@cached_property
	def components(self):
		"""
		Extends components with the external moment.
		"""

		components = {}

		# every parameter and coordinate is passed explicitly - otherwise the
		# subsystem would build its own symbols and the external moment could
		# act on coordinates that do not belong to the rotor
		self._elastic_rotor = FreeEngineMDOF(
			m=self.m,
			I_engine=self.I_engine,
			e=self.e,
			k=self.k,
			h=self.h,
			v=self.v,
			phi=self.phi,
			ivar=self.ivar,
		)(label="Rotor with Shaft")

		self._external_moment_comp = Force(
			self.M_t,  # Value of the moment
			pos1=self.phi,  # The coordinate on which the moment acts
			qs=self.qs,
		)(label="External moment")

		components["_elastic_rotor"] = self._elastic_rotor
		components["_external_moment"] = self._external_moment_comp

		return components

	def symbols_description(self):
		"""
		Descriptions of every symbol appearing in the equations of motion,
		including the derivatives - the report tables are built from the
		equations, so the derivatives need entries of their own.
		"""

		t = self.ivar

		self.sym_desc_dict = {
			self.m: r"Mass of the rotor",
			self.I_engine: r"Mass moment of inertia of the rotor",
			self.e: r"Eccentricity of the centre of mass",
			self.k: r"Transverse stiffness of the shaft",
			self.M_t: r"External disturbance torque",
			self.h: r"Horizontal displacement of the mass centre",
			self.h.diff(t): r"Horizontal velocity of the mass centre",
			self.h.diff(t, 2): r"Horizontal acceleration of the mass centre",
			self.v: r"Vertical displacement of the mass centre",
			self.v.diff(t): r"Vertical velocity of the mass centre",
			self.v.diff(t, 2): r"Vertical acceleration of the mass centre",
			self.phi: r"Angular position of the rotor",
			self.phi.diff(t): r"Angular velocity of the rotor",
			self.phi.diff(t, 2): r"Angular acceleration of the rotor",
			t: r"Time",
		}

		return self.sym_desc_dict

	def unit_dict(self):

		t = self.ivar

		units_dict = {
			self.m: ureg.kilogram,
			self.I_engine: ureg.kilogram * ureg.meter**2,
			self.e: ureg.meter,
			self.k: ureg.newton / ureg.meter,
			self.M_t: ureg.newton * ureg.meter,
			self.h: ureg.meter,
			self.h.diff(t): ureg.meter / ureg.second,
			self.h.diff(t, 2): ureg.meter / ureg.second**2,
			self.v: ureg.meter,
			self.v.diff(t): ureg.meter / ureg.second,
			self.v.diff(t, 2): ureg.meter / ureg.second**2,
			self.phi: ureg.radian,
			self.phi.diff(t): ureg.radian / ureg.second,
			self.phi.diff(t, 2): ureg.radian / ureg.second**2,
			t: ureg.second,
		}

		return units_dict
