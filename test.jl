include("simulate.jl")
clearconsole()
I1 = [1 0 0; 0 0.1 0; 0 0 1.2]
m = 1.0
using Revise
# Test the code.
# for i in 1:100
#     #ang0 = [45*pi/180;60*pi/180;30*pi/180];
#     ang0 = rand(3);
#     q = angle2quat(ang0,"321");
#     ang = quat2angle(q,"321");
#     e = ang - ang0;
#     println("Norm error:", norm(e));
# end

# Simulate rigid body dynamics
# We have to now call ode solver with parameters
# -- this will have a function that will be called to compute all time varying
# -- quantities.


# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,"quaternions")
RbI = InertialFrameAsRB()
x0R1 = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
R1 = initialiseRigidBody(R1,x0R1)
axis = [0.0 1.0 0.0][:]
rj1 = [0.0 0.0 0.0][:]
rj2 = [0.0 0.0 0.0][:]
# x0R2 = rand(14)
# R2 = RigidBody(m,I1,"quaternions")
# R2 = initialiseRigidBody(R1,x0R2)
j = Joint(RbI,R1,rj1,rj2,"Revolute",axis)
tEnd = 1.0
tSpan = 0.01
simulate(tEnd,tSpan,j)
