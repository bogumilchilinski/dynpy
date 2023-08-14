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
from .trolley import SpringMassSystem, ForcedSpringMassSystem, SpringDamperMassSystem, SprungTrolley, NonLinearTrolley, TrolleyWithPendulum, DampedTrolleyWithPendulum, ForcedNonLinearTrolley, ForcedTrolleyWithSpring, ForcedDampedTrolleyWithSpring, ForcedDampedDisconnectedTrolleysWithSprings, VariableMassTrolleyWithPendulum, VariableMassTrolleyWithPendulumRayleighDamping, TrolleyWithTMD, VariableMassTrolleyWithPendulumFunction
from .trolleys import DDoFTwoNonLinearTrolleys, TrolleysWithSprings, ForcedTrolleysWithSprings, ForcedDampedTrolleysWithSprings, ForcedThreeTrolleysWithSprings, ForcedTrolleysWithNonLinearSprings, ForcedDisconnectedTrolleysWithSprings, ForcedDampedThreeTrolleysWithSprings
from .vehicle import UndampedVehicleSuspension, DampedVehicleSuspension, UndampedSymmetricalVehicleSuspension, QuarterOfVehicle
from .shaft import SDoFShaft, SDoFShaftHarmonicExcitation, DoubleDiskShaft, DoubleDiskShaftHarmonicExcitation,  DampedDoubleDiskShaftHarmonicExcitation, SDoFDampedShaft
from .pendulum import Pendulum, PulledPendulum, FreePendulum, KinematicallyExcitedInvertedPendulum, MDoFElasticPendulum, ForcedTriplePendulum, Winch, DDOFCoupledPendulum, DDOFLinearizedCoupledPendulum, InvertedPendulumDDoF, DampedElasticPendulum, LinearizedTriplePendulum

from .skyscraper import skyscraper
from .engine import EngineHousing, FreeEngine,  FreeEngineWithConstrains, FreeEngineDDOF, Engine, EngineVerticalSpringGravity,  EngineConstantVelocityVerticalSpringGravity,  EngineConstantWithReactionForce,  DampedEngineVerticalSpringGravity,  DampedEngineConstantVelocityVerticalSpringGravity,  InlineEnginePerpendicularSprings,  NonLinearInlineEnginePerpendicularSpringsGravity,  NonLinearBoxerEnginePerpendicularSprings,  DampedEngine, NonlinearEngine,  StraightNonlinearEngine, EngineWithTMD,  VeeEnginePerpendicularSprings
from .disk import ForcedNonLinearDisc, ForcedNonLinearDiscSpring
from .disks import TwoForcedNonLinearDisks, TwoDisksWithThreeSprings, ForcedDisksWithSerialSprings, ForcedDisksWithParallelSprings, MDoFForcedSimpleDisksWithSerialSprings, MDoFForcedSimpleDisksWithParallelSprings
from .tensioner import BlowerToothedBelt, DampedBlowerToothedBelt
from .tmd import TunedMassDamper, TunedMassDamperRelativeMotion

from .projectile import MissileTrajectoryAirless, MissileTrajectory, MissileTrajectoryAerodynamic
from .tmac import SDOFWinchSystem, SDOFWinchSystemTest, SDOFDrivetrainVehicleSystem

from .principles import DampedMeasuringTool, KinematicClutchWithSprings, LagrangeIBlocksOnInclinedPlane, MaterialPointMovement

from .offshore import CompoundSystem, SDOFSemiSubmersiblePlatform, DDOFSemiSubmersiblePlatform

#from .car import CarMovementPIDAdjust, GearboxEngine, CarMovementRegulatedThrottleACC, CarMovementRegulatedThrottle, CarMovementAdjustableThrottle, OwnCombustionEngine, SecondGearCarMovement, ThirdGearCarMovement, FourthGearCarMovement, FifthGearCarMovement, CarMovementConstantThrottle

#from .motorcoach import FourDOFTrolleySuspension, DDOFTrolleySuspension2
#TODO
#from .motorcoach import DDOFTrolleySuspension

#from .chair import DampedChair4DOF, DampedChairDDOF, DampedChairSimplifiedDDOF2,


#TODO
#from .disks import ForcedDisksWithParallelSprings2, MDoFForcedDisksWithParallelSprings
#from .engine import BoxerEnginePerpendicularSprings, NonLinearVeeEnginePerpendicularSprings
#from .motorcoach import FourDOFTrolleySuspension, DDOFTrolleySuspension, DDOFTrolleySuspension2
#from .pendulum import ExcitedPendulum, DampedPendulum, ExcitedDampedPendulum, MDoFLinearizedThreePendulumsWithSprings, , ThreePendulumsWithSprings, DDoFWinch
#from .principles import  LagrangeIOnMathFunction
#from .shaft import TripleShaft

#from .trolley import DoubleTrolleyWithNonlinearSprings, DoubleTrolleyDifferentWheels, MDoFDoubleTrolleyWithNonlinearSprings, SDoFDoubleTrolleyDifferentWheels, MDoFDoubleTrolleyDifferentWheels, MDoFTMD 
#from .trolleys import DampedTrolleysWithSprings, DampedTrolleysWithSprings, ForcedTrolleysWithSprings, MDoFForcedTrolleysWithSprings, MDoFDampedTrolleysWithSprings
#from .vehicle import DampedSymmetricalVehicleSuspension, SimplifySuspension




