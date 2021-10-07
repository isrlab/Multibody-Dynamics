# include("../src/plotSol.jl")
# include("../src/OrientationConversion.jl")
# include("../src/simulate.jl")
#
# using JLD
# clearconsole()
#
# savedSol = load("ModRob1dSecSimRK45NoTol.jld") # DifferentialEquations
# # sol = load("ModRob1dSecRK4.jld") # Native RK
# solFinal = savedSol["solFinal"]
# tSim = savedSol["tSim"]

## Assigning solutions
solBody1 = solFinal[1]
solGimbal1 = solFinal[2]
solProp11 = solFinal[3]
solProp12 = solFinal[4]
solBox = solFinal[5]
solBody2 = solFinal[6]
solGimbal2 = solFinal[7]
solProp21 = solFinal[8]
solProp22 = solFinal[9]

## Plotting
jointLoc2 = Matrix{Float64}(undef,length(tSim),3)
jointLoc3 = Matrix{Float64}(undef,length(tSim),3)
jointLoc4 = Matrix{Float64}(undef,length(tSim),3)
jointLoc5 = Matrix{Float64}(undef,length(tSim),3)
jointLoc6 = Matrix{Float64}(undef,length(tSim),3)
jointLoc7 = Matrix{Float64}(undef,length(tSim),3)
jointLoc8 = Matrix{Float64}(undef,length(tSim),3)
jointLoc9 = Matrix{Float64}(undef,length(tSim),3)

ωBody1 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal1 = Matrix{Float64}(undef,length(tSim),3)
ωProp11 = Matrix{Float64}(undef,length(tSim),3)
ωProp12 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal1InBody1 = Matrix{Float64}(undef,length(tSim),3)
ωRel1 = Matrix{Float64}(undef,length(tSim),3)
ωRel2 = Matrix{Float64}(undef,length(tSim),3)
ωRel3 = Matrix{Float64}(undef,length(tSim),3)

ωBox = Matrix{Float64}(undef,length(tSim),3)

ωBody2 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal2 = Matrix{Float64}(undef,length(tSim),3)
ωProp21 = Matrix{Float64}(undef,length(tSim),3)
ωProp22 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal2InBody2 = Matrix{Float64}(undef,length(tSim),3)
ωRel4 = Matrix{Float64}(undef,length(tSim),3)
ωRel5 = Matrix{Float64}(undef,length(tSim),3)
ωRel6 = Matrix{Float64}(undef,length(tSim),3)
for i=1:length(tSim)
    Body1.dcm = quat2dcm(solBody1.β[i,:])
    Gimbal1.dcm = quat2dcm(solGimbal1.β[i,:])
    Prop11.dcm = quat2dcm(solProp11.β[i,:])
    Prop12.dcm = quat2dcm(solProp12.β[i,:])
    Box.dcm = quat2dcm(solBox.β[i,:])

    ωBox[i,:] = angVel(solBox.β[i,:],solBox.βdot[i,:])
    ωBody1[i,:] = angVel(solBody1.β[i,:],solBody1.βdot[i,:])
    ωGimbal1[i,:] = angVel(solGimbal1.β[i,:],solGimbal1.βdot[i,:])
    ωProp11[i,:] = angVel(solProp11.β[i,:],solProp11.βdot[i,:])
    ωProp12[i,:] = angVel(solProp12.β[i,:],solProp12.βdot[i,:])

    ωGimbal1InBody1[i,:] = quat2dcm(solBody1.β[i,:])*
                         transpose(quat2dcm(solGimbal1.β[i,:]))*
                         ωGimbal1[i,:]
    ωRel1[i,:] = ωGimbal1InBody1[i,:] - ωBody1[i,:]
    ωRel2[i,:] = quat2dcm(solGimbal1.β[i,:])*
                 transpose(quat2dcm(solProp11.β[i,:]))*
                 ωProp11[i,:] - ωGimbal1[i,:]
    ωRel3[i,:] = quat2dcm(solGimbal1.β[i,:])*
                 transpose(quat2dcm(solProp12.β[i,:]))*
                 ωProp12[i,:] - ωGimbal1[i,:]

    jointLoc2[i,:] = solBody1.r[i,:] + transpose(Body1.dcm)*rjBody  -                    solGimbal1.r[i,:] -                     transpose(Gimbal1.dcm)*rjGimbal1
    jointLoc3[i,:] = solGimbal1.r[i,:] + transpose(Gimbal1.dcm)*rjGimbal2 - solProp11.r[i,:] - transpose(Prop11.dcm)*rjProp1
    jointLoc4[i,:] = solGimbal1.r[i,:] +
                     transpose(Gimbal1.dcm)*rjGimbal3 -
                     solProp12.r[i,:] - transpose(Prop12.dcm)*rjProp2
    jointLoc5[i,:] = solBody1.r[i,:] + transpose(Body1.dcm)*rjCac1 - solBox.r[i,:] - transpose(Box.dcm)*rjBox1
end
for i=1:length(tSim)
     Body2.dcm = quat2dcm(solBody2.β[i,:])
     Gimbal2.dcm = quat2dcm(solGimbal2.β[i,:])
     Prop21.dcm = quat2dcm(solProp21.β[i,:])
     Prop22.dcm = quat2dcm(solProp22.β[i,:])
     Box.dcm = quat2dcm(solBox.β[i,:])

     ωBody2[i,:] = angVel(solBody2.β[i,:],solBody2.βdot[i,:])
     ωGimbal2[i,:] = angVel(solGimbal2.β[i,:],solGimbal2.βdot[i,:])
     ωProp21[i,:] = angVel(solProp21.β[i,:],solProp21.βdot[i,:])
     ωProp22[i,:] = angVel(solProp22.β[i,:],solProp22.βdot[i,:])

     ωGimbal2InBody2[i,:] = quat2dcm(solBody2.β[i,:])*
                          transpose(quat2dcm(solGimbal2.β[i,:]))*
                          ωGimbal2[i,:]
     ωRel4[i,:] = ωGimbal2InBody2[i,:] - ωBody2[i,:]
     ωRel5[i,:] = quat2dcm(solGimbal2.β[i,:])*
                  transpose(quat2dcm(solProp21.β[i,:]))*
                  ωProp21[i,:] - ωGimbal2[i,:]
     ωRel6[i,:] = quat2dcm(solGimbal2.β[i,:])*
                  transpose(quat2dcm(solProp22.β[i,:]))*
                  ωProp22[i,:] - ωGimbal2[i,:]

     jointLoc7[i,:] = solBody2.r[i,:] + transpose(Body2.dcm)*rjBody  -                    solGimbal2.r[i,:] -                     transpose(Gimbal2.dcm)*rjGimbal1
     jointLoc8[i,:] = solGimbal2.r[i,:] + transpose(Gimbal2.dcm)*rjGimbal2 - solProp21.r[i,:] - transpose(Prop21.dcm)*rjProp1
     jointLoc9[i,:] = solGimbal2.r[i,:] +
                      transpose(Gimbal2.dcm)*rjGimbal3 -
                      solProp22.r[i,:] - transpose(Prop22.dcm)*rjProp2
     jointLoc6[i,:] = solBody2.r[i,:] + transpose(Body2.dcm)*rjCac2 - solBox.r[i,:] - transpose(Box.dcm)*rjBox2
end

## Plotting
# Joint Locations
plotPos(tSim,jointLoc2)
plotPos(tSim,jointLoc3)
plotPos(tSim,jointLoc4)
plotPos(tSim,jointLoc5)
plotPos(tSim,jointLoc6)
plotPos(tSim,jointLoc7)
plotPos(tSim,jointLoc8)
plotPos(tSim,jointLoc9)

plotErrNorm(tSim,solBody1.β)
plotErrNorm(tSim,solGimbal1.β)
plotErrNorm(tSim,solProp11.β)
plotErrNorm(tSim,solProp12.β)
plotErrNorm(tSim,solBody2.β)
plotErrNorm(tSim,solGimbal2.β)
plotErrNorm(tSim,solProp21.β)
plotErrNorm(tSim,solProp22.β)

plotAngVel(tSim,ωBody1)
plotAngVel(tSim,ωGimbal1)
plotAngVel(tSim,ωProp11)
plotAngVel(tSim,ωProp12)
plotAngVel(tSim,ωRel1)
plotAngVel(tSim,ωRel2)
plotAngVel(tSim,ωRel3)

plotAngVel(tSim,ωBox)

plotAngVel(tSim,ωBody2)
plotAngVel(tSim,ωGimbal2)
plotAngVel(tSim,ωProp21)
plotAngVel(tSim,ωProp22)
plotAngVel(tSim,ωRel4)
plotAngVel(tSim,ωRel5)
plotAngVel(tSim,ωRel6)
