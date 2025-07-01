"""
This module provides the examples of mechanical models being practically used in calculations
"""

import importlib

from .bridge import (
    BeamBridge,
    BeamBridgeDamped,
    BeamBridgeDampedTMD,
    BeamBridgeTMD,
    BeamElasticity,
    BeamStructure,
)
from .car import (
    CarMovementAdjustableThrottle,
    CarMovementConstantThrottle,
    CarMovementPIDAdjust,
    CarMovementRegulatedThrottle,
    CarMovementRegulatedThrottleACC,
    FifthGearCarMovement,
    FourthGearCarMovement,
    GearboxEngine,
    OwnCombustionEngine,
    SecondGearCarMovement,
    ThirdGearCarMovement,
)
from .disk import (
    ForcedNonLinearDisk,
    ForcedNonLinearDiskSpring,
    ForcedRollingHalfDisk,
    RollingHalfDisk,
)
from .disks import (
    ForcedDisksWithParallelSprings,
    ForcedDisksWithSerialSprings,
    MDoFForcedDisksWithParallelSprings,
    MDoFForcedSimpleDisksWithParallelSprings,
    MDoFForcedSimpleDisksWithSerialSprings,
    TwoDisksWithThreeSprings,
    TwoForcedNonLinearDisks,
)
from .engine import (
    BoxerEnginePerpendicularSprings,
    DampedEngine,
    DampedEngineConstantVelocityVerticalSpringGravity,
    DampedEngineVerticalSpringGravity,
    Engine,
    EngineConstantVelocityVerticalSpringGravity,
    EngineConstantWithReactionForce,
    EngineHousing,
    EngineVerticalSpringGravity,
    EngineWithTMD,
    FreeEngine,
    FreeEngineDDOF,
    FreeEngineWithConstrains,
    InlineEnginePerpendicularSprings,
    NonLinearBoxerEnginePerpendicularSprings,
    NonlinearEngine,
    NonLinearInlineEnginePerpendicularSpringsGravity,
    NonLinearVeeEnginePerpendicularSprings,
    StraightNonlinearEngine,
    VeeEnginePerpendicularSprings,
)
from .offshore import (
    CompoundSystem,
    DDOFSemiSubmersiblePlatform,
    SDOFSemiSubmersiblePlatform,
)
from .pendulum import (
    CompoundPendulum,
    DampedElasticPendulum,
    DampedMeasuringTool,
    DampedPendulum,
    DDOFCoupledPendulum,
    DDOFLinearizedCoupledPendulum,
    DDoFWinch,
    DoublePendulum,
    ExcitedDampedPendulum,
    ExcitedPendulum,
    ForcedCompoundPendulum,
    ForcedRollingBar,
    ForcedTriplePendulum,
    FreePendulum,
    InvertedPendulumDDoF,
    KinematicallyExcitedInvertedPendulum,
    LinearizedTriplePendulum,
    MDoFElasticPendulum,
    MDoFLinearizedThreePendulumsWithSprings,
    Pendulum,
    PendulumKinematicExct,
    PulledPendulum,
    RollingBar,
    Winch,
)
from .principles import *
from .principles import (
    KinematicClutchWithSprings,
    LagrangeIBlocksOnInclinedPlane,
    MaterialPointMovement,
)
from .projectile import (
    MissileTrajectory,
    MissileTrajectoryAerodynamic,
    MissileTrajectoryAirless,
)
from .shaft import (
    DampedDoubleDiskShaftHarmonicExcitation,
    DoubleDiskShaft,
    DoubleDiskShaftHarmonicExcitation,
    SDoFDampedShaft,
    SDoFShaft,
    SDoFShaftHarmonicExcitation,
)
from .skyscraper import skyscraper
from .tensioner import BlowerToothedBelt, DampedBlowerToothedBelt
from .tmac import SDOFDrivetrainVehicleSystem, SDOFWinchSystem, SDOFWinchSystemTest
from .tmd import TunedMassDamper, TunedMassDamperRelativeMotion
from .trolley import (  # , ForcedDampedTrolleyWithSpring
    DampedTrolleyWithPendulum,
    ForcedNonLinearTrolley,
    ForcedSpringMassSystem,
    ForcedTrolleyWithSpring,
    NonLinearTrolley,
    SpringDamperMassSystem,
    SpringMassSystem,
    SprungTrolley,
    TrolleyWithPendulum,
    TrolleyWithTMD,
    VariableMassTrolleyWithPendulum,
    VariableMassTrolleyWithPendulumFunction,
    VariableMassTrolleyWithPendulumRayleighDamping,
)
from .trolleys import (
    DDoFTwoNonLinearTrolleys,
    ForcedDampedDisconnectedTrolleysWithSprings,
    ForcedDampedThreeTrolleysWithSprings,
    ForcedDampedTrolleysWithSprings,
    ForcedDisconnectedTrolleysWithSprings,
    ForcedThreeTrolleysWithSprings,
    ForcedTrolleysWithNonLinearSprings,
    ForcedTrolleysWithSprings,
    TrolleysWithSprings,
)
from .vehicle import (
    DampedVehicleSuspension,
    QuarterOfVehicle,
    UndampedSymmetricalVehicleSuspension,
    UndampedVehicleSuspension,
)

# from .bridge import *
# from .principles import *
# from .trolley import *
# from .vehicle import *
# from .shaft import *
# from .pendulum import *
# from .car import *
# from .skyscraper import *
# from .engine import *
# from .disk import *
# from .tensioner import *
# from .tmd import *
# from .chair import *
# from .projectile import *
# from .tmac import *
# from .motorcoach import *
# from .offshore import *


# from .motorcoach import FourDOFTrolleySuspension, DDOFTrolleySuspension2
# TODO
# from .motorcoach import DDOFTrolleySuspension

# from .chair import DampedChair4DOF, DampedChairDDOF, DampedChairSimplifiedDDOF2,


# TODO
# from .engine import #163,#162 (TODO)
# from .motorcoach import FourDOFTrolleySuspension, DDOFTrolleySuspension, DDOFTrolleySuspension2
# from .principles import  LagrangeIOnMathFunction
# from .shaft import TripleShaft
# from .trolleys import DampedTrolleysWithSprings, DampedTrolleysWithSprings, ForcedTrolleysWithSprings, MDoFForcedTrolleysWithSprings, MDoFDampedTrolleysWithSprings, DoubleTrolleyWithNonlinearSprings, DoubleTrolleyDifferentWheels, MDoFDoubleTrolleyWithNonlinearSprings, SDoFDoubleTrolleyDifferentWheels, MDoFDoubleTrolleyDifferentWheels, MDoFTMD
# from .vehicle import DampedSymmetricalVehicleSuspension, SimplifySuspension
