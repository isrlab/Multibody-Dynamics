# Attempting to use finitediff instead of autodiff for linearizations
using FiniteDifferences
x0, u0 = getXU_0(j);

m = central_fdm(10,1)
A1 = FiniteDifferences.jacobian(m, z->fxdot(z,u0,j,g),x0)[1]
B1 = FiniteDifferences.jacobian(m, z->fxdot(x0,z,j,g),u0)[1]

A2, B2 = linearizeDiff(x0,u0,j)
##
println()
println("Aerr = ", norm(A1-A2))

println("Berr = ", norm(B1 - B2))
