# function that returns the joint tree for an n-rotor vtol vehicle
include("simulate.jl")
function gen_nRotor(n::Integer)
    # Inertial and mass properties from:
    # Hossain, M. Raju, and Nicholas Krouglicof. "Multi-body dynamics modeling & control of quadrotor helicopter using bond graph." Proceedings of the International Conference on Bond Graph Modeling and Simulation. 2016.

    mb = 0.30155; mp = 0.011721;

    Ib = zeros(3,3);
    Ib[1,1] = 6.4e-4; Ib[2,2] = 6.4e-4; Ib[3,3] = 1.22e-3;

    Ip = zeros(3,3);
    Ip[1,1] = 6.4061e-7; Ip[2,2] = 6.2732e-6; Ip[3,3] = 3.008e-5;

    InFrame = InertialFrameAsRB()
    QuadFrame = RigidBody(mb,Ib,2)
    x0Quad = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)];
    initialiseRigidBody!(QuadFrame,x0Quad);

    Props = Vector{RigidBody}(undef,n);
    j = Vector{Joint}(undef, n+1);
    j[1] = Joint(InFrame,QuadFrame,zeros(3),zeros(3));

    th = 2*pi/n;
    l = 1.0;
    axis = [0.0 0.0 1.0][:];
    for i=1:n
        Props[i] = RigidBody(mp, Ip, 2+i);
        ang =(i-1)*th; # angle
        r0Prop_i = [l*cos(ang); l*sin(ang); 0.0];
        β0Prop_i = abs.([cos(ang/2); sin(ang/2).*[0.0;0.0;1.0]]);
        x0Prop_i = [r0Prop_i; β0Prop_i; zeros(7)];
        initialiseRigidBody!(Props[i],x0Prop_i);


        j[i+1] = Joint(QuadFrame, Props[i], r0Prop_i, zeros(3), type = "Revolute", axis = axis);
    end
    return j
end
## simulate quadcopter
# clearconsole();
# qRotor_j = gen_nRotor(4);
#
# tEnd = 0.1;
# tSpan = 0.01;
# g =([0.0,0.0,-9.81]); # Gravity Vector.
# tSim, solFinal = simulate(tEnd, tSpan, qRotor_j, g=g);
#
# solQuad = solFinal[1];
# solProp1 = solFinal[2];
# solProp2 = solFinal[3];
# solProp3 = solFinal[4];
# solProp4 = solFinal[5];
