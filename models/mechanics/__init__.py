"""
This module provides the examples of mechanical models being practically used in calculations
"""

import importlib

#from .bridge import *
#from .principles import *
#from .trolley import *
#from .vehicle import *
#from .shaft import *
#from .pendulum import *
#from .car import *
#from .skyscraper import *
#from .engine import *
#from .disk import *
#from .tensioner import *
#from .tmd import *
#from .chair import *
#from .projectile import *
#from .tmac import *
#from .motorcoach import *
#from .offshore import *

from dynpy.models.mechanics.bridge import BeamStructure, BeamElasticity, BeamBridge, BeamBridgeDamped, BeamBridgeTMD,  BeamBridgeDampedTMD
from dynpy.models.mechanics.principles import *
from dynpy.models.mechanics.trolley import SpringMassSystem, ForcedSpringMassSystem, SpringDamperMassSystem, SprungTrolley, NonLinearTrolley, TrolleyWithPendulum, DampedTrolleyWithPendulum, ForcedNonLinearTrolley, TrolleysWithSprings, ForcedTrolleysWithSprings,  ForcedDampedTrolleysWithSprings, ForcedTrolleysWithNonLinearSprings, ForcedTrolleyWithSpring, ForcedDampedTrolleyWithSpring, ForcedThreeTrolleysWithSprings, ForcedDisconnectedTrolleysWithSprings, ForcedDampedThreeTrolleysWithSprings, ForcedDampedDisconnectedTrolleysWithSprings, DDoFTwoNonLinearTrolleys, VariableMassTrolleyWithPendulum, VariableMassTrolleyWithPendulumRayleighDamping, TrolleyWithTMD, VariableMassTrolleyWithPendulumFunction
from dynpy.models.mechanics.vehicle import UndampedVehicleSuspension, DampedVehicleSuspension, UndampedSymmetricalVehicleSuspension, QuarterOfVehicle
from dynpy.models.mechanics.shaft import SDoFShaft, SDoFShaftHarmonicExcitation, DoubleDiskShaft, DoubleDiskShaftHarmonicExcitation,  DampedDoubleDiskShaftHarmonicExcitation, SDoFDampedShaft
from dynpy.models.mechanics.pendulum import Pendulum, PulledPendulum, FreePendulum, KinematicallyExcitedInvertedPendulum, MDoFElasticPendulum, ForcedTriplePendulum, Winch, DDOFCoupledPendulum, DDOFLinearizedCoupledPendulum, InvertedPendulumDDoF, DampedElasticPendulum
from dynpy.models.mechanics.car import *
from dynpy.models.mechanics.skyscraper import skyscraper
from dynpy.models.mechanics.engine import EngineHousing, FreeEngine,  FreeEngineWithConstrains, FreeEngineDDOF, Engine, EngineVerticalSpringGravity,  EngineConstantVelocityVerticalSpringGravity,  EngineConstantWithReactionForce,  DampedEngineVerticalSpringGravity,  DampedEngineConstantVelocityVerticalSpringGravity,  InlineEnginePerpendicularSprings,  NonLinearInlineEnginePerpendicularSpringsGravity,  NonLinearBoxerEnginePerpendicularSprings,  DampedEngine, NonlinearEngine,  StraightNonlinearEngine, EngineWithTMD,  VeeEnginePerpendicularSprings
from dynpy.models.mechanics.disk import ForcedNonLinearDisc, ForcedNonLinearDiscSpring, TwoForcedNonLinearDisks, TwoDisksWithThreeSprings
from dynpy.models.mechanics.tensioner import BlowerToothedBelt, DampedBlowerToothedBelt
from dynpy.models.mechanics.tmd import TunedMassDamper, TunedMassDamperRelativeMotion
from dynpy.models.mechanics.chair import *
from dynpy.models.mechanics.projectile import MissileTrajectoryAirless, MissileTrajectory, MissileTrajectoryAerodynamic
from dynpy.models.mechanics.tmac import SDOFWinchSystem, SDOFWinchSystemTest
from dynpy.models.mechanics.motorcoach import *
from dynpy.models.mechanics.principles import DampedMeasuringTool, KinematicClutchWithSprings
from dynpy.models.mechanics.offshore import CompoundSystem, SDOFSemiSubmersiblePlatform, DDOFSemiSubmersiblePlatform


#TODO
#from dynpy.models.mechanics.disk import ForcedDisksWithSerialSprings, ForcedDisksWithParallelSprings, ForcedDisksWithParallelSprings2, MDoFForcedSimpleDisksWithSerialSprings, MDoFForcedSimpleDisksWithParallelSprings, MDoFForcedDisksWithParallelSprings
#from dynpy.models.mechanics.engine import BoxerEnginePerpendicularSprings, NonLinearVeeEnginePerpendicularSprings
#from dynpy.models.mechanics.motorcoach import FourDOFTrolleySuspension, DDOFTrolleySuspension, DDOFTrolleySuspension2
#from dynpy.models.mechanics.pendulum import ExcitedPendulum, DampedPendulum, ExcitedDampedPendulum, MDoFLinearizedThreePendulumsWithSprings, LinearizedTriplePendulum, ThreePendulumsWithSprings, DDoFWinch
#from dynpy.models.mechanics.principles import LagrangeIBlocksOnInclinedPlane, LagrangeIOnMathFunction, MaterialPointMovement
#from dynpy.models.mechanics.shaft import TripleShaft
#from dynpy.models.mechanics.tmac import SDOFDrivetrainVehicleSystem
#from dynpy.models.mechanics.trolley import DampedTrolleysWithSprings, DoubleTrolleyWithNonlinearSprings, DampedTrolleysWithSprings, DoubleTrolleyDifferentWheels, ForcedTrolleysWithSprings, MDoFDoubleTrolleyWithNonlinearSprings, MDoFDampedTrolleysWithSprings, SDoFDoubleTrolleyDifferentWheels, MDoFDoubleTrolleyDifferentWheels, MDoFForcedTrolleysWithSprings, MDoFTMD 
#from dynpy.models.mechanics.vehicle import DampedSymmetricalVehicleSuspension, SimplifySuspension




