module UdKaDyn

include("RigidBody.jl")
include("Joint.jl")
include("OrientationConversion.jl")
include("extF.jl")
include("simulate.jl")
include("Force.jl")

include("linearize.jl")
include("nRotor.jl")
include("plotSol.jl")
# include("PointMass.jl")
include("trim_FinDiff.jl")
include("trim_kron.jl")
include("trim_kronLazy.jl")

export gen_nRotor, simulate, RigidBody, InertialFrameAsRB, Joint, checkRevJointIn, plotErrNorm, plotPos, plotQuat, plotVel, plotAngVel, linearize, initialiseRigidBody!
end # module
