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

from .bridge import BeamStructure, BeamElasticity, BeamBridge, BeamBridgeDamped, BeamBridgeTMD,  BeamBridgeDampedTMD
from .principles import *
from .trolley import SpringMassSystem, ForcedSpringMassSystem, SpringDamperMassSystem, SprungTrolley, NonLinearTrolley, TrolleyWithPendulum, DampedTrolleyWithPendulum, ForcedNonLinearTrolley, TrolleysWithSprings, ForcedTrolleysWithSprings,  ForcedDampedTrolleysWithSprings, ForcedTrolleysWithNonLinearSprings, ForcedTrolleyWithSpring, ForcedDampedTrolleyWithSpring, ForcedThreeTrolleysWithSprings, ForcedDisconnectedTrolleysWithSprings, ForcedDampedThreeTrolleysWithSprings, ForcedDampedDisconnectedTrolleysWithSprings, DDoFTwoNonLinearTrolleys, VariableMassTrolleyWithPendulum, VariableMassTrolleyWithPendulumRayleighDamping, TrolleyWithTMD, VariableMassTrolleyWithPendulumFunction
from .vehicle import UndampedVehicleSuspension, DampedVehicleSuspension, UndampedSymmetricalVehicleSuspension, QuarterOfVehicle
from .shaft import SDoFShaft, SDoFShaftHarmonicExcitation, DoubleDiskShaft, DoubleDiskShaftHarmonicExcitation,  DampedDoubleDiskShaftHarmonicExcitation, SDoFDampedShaft
from .pendulum import Pendulum, PulledPendulum, FreePendulum, KinematicallyExcitedInvertedPendulum, MDoFElasticPendulum, ForcedTriplePendulum, Winch, DDOFCoupledPendulum, DDOFLinearizedCoupledPendulum, InvertedPendulumDDoF, DampedElasticPendulum
from .car import *
from .skyscraper import skyscraper
from .engine import EngineHousing, FreeEngine,  FreeEngineWithConstrains, FreeEngineDDOF, Engine, EngineVerticalSpringGravity,  EngineConstantVelocityVerticalSpringGravity,  EngineConstantWithReactionForce,  DampedEngineVerticalSpringGravity,  DampedEngineConstantVelocityVerticalSpringGravity,  InlineEnginePerpendicularSprings,  NonLinearInlineEnginePerpendicularSpringsGravity,  NonLinearBoxerEnginePerpendicularSprings,  DampedEngine, NonlinearEngine,  StraightNonlinearEngine, EngineWithTMD,  VeeEnginePerpendicularSprings
from .disk import ForcedNonLinearDisc, ForcedNonLinearDiscSpring, TwoForcedNonLinearDisks, TwoDisksWithThreeSprings
from .tensioner import BlowerToothedBelt, DampedBlowerToothedBelt
from .tmd import TunedMassDamper, TunedMassDamperRelativeMotion

from .projectile import MissileTrajectoryAirless, MissileTrajectory, MissileTrajectoryAerodynamic
from .tmac import SDOFWinchSystem, SDOFWinchSystemTest

from .principles import DampedMeasuringTool, KinematicClutchWithSprings

from .offshore import CompoundSystem, SDOFSemiSubmersiblePlatform, DDOFSemiSubmersiblePlatform
from .motorcoach import *
from .chair import *


#TODO
#from .disk import ForcedDisksWithSerialSprings, ForcedDisksWithParallelSprings, ForcedDisksWithParallelSprings2, MDoFForcedSimpleDisksWithSerialSprings, MDoFForcedSimpleDisksWithParallelSprings, MDoFForcedDisksWithParallelSprings
#from .engine import BoxerEnginePerpendicularSprings, NonLinearVeeEnginePerpendicularSprings
#from .motorcoach import FourDOFTrolleySuspension, DDOFTrolleySuspension, DDOFTrolleySuspension2
#from .pendulum import ExcitedPendulum, DampedPendulum, ExcitedDampedPendulum, MDoFLinearizedThreePendulumsWithSprings, LinearizedTriplePendulum, ThreePendulumsWithSprings, DDoFWinch
#from .principles import LagrangeIBlocksOnInclinedPlane, LagrangeIOnMathFunction, MaterialPointMovement
#from .shaft import TripleShaft
#from .tmac import SDOFDrivetrainVehicleSystem
#from .trolley import DampedTrolleysWithSprings, DoubleTrolleyWithNonlinearSprings, DampedTrolleysWithSprings, DoubleTrolleyDifferentWheels, ForcedTrolleysWithSprings, MDoFDoubleTrolleyWithNonlinearSprings, MDoFDampedTrolleysWithSprings, SDoFDoubleTrolleyDifferentWheels, MDoFDoubleTrolleyDifferentWheels, MDoFForcedTrolleysWithSprings, MDoFTMD 
#from .vehicle import DampedSymmetricalVehicleSuspension, SimplifySuspension




