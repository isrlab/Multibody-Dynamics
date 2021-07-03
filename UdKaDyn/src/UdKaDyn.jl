module UdKaDyn

include("extF.jl")
include("Force.jl")
include("Joint.jl")
include("linearize.jl")
include("nRotor.jl")
include("OrientationConversion.jl")
include("plotSol.jl")
include("PointMass.jl")
include("simulate.jl")
include("trim_FinDiff.jl")
include("trim_kron.jl")
include("trim_kronLazy.jl")

export gen_nRotor, simulate, RigidBody, InertialFrameAsRB, Joint, checkRevJointIn, plotErrNorm, plotPos, plotQuat, plotVel, plotAngVel, linearize
end # module
