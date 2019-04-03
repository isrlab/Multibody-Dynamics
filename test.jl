include("RigidBody.jl");
I1 = [1 0 0; 0 0.1 0; 0 0 1.2];
m = 1.0;

R = RigidBody(m,I1,"321");

include("OrientationConversion.jl");

# Test the code.
for i in 1:100
    #ang0 = [45*pi/180;60*pi/180;30*pi/180];
    ang0 = rand(3);
    q = angle2quat(ang0,"321");
    ang = quat2angle(q,"321");
    e = ang - ang0;
    println("Norm error:", norm(e));
end

# Simulate rigid body dynamics
# We have to now call ode solver with parameters
# -- this will have a function that will be called to compute all time varying
# -- quantities.
