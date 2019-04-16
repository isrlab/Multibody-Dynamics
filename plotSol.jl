using Plots
pyplot();

function plotPos(tSim,rSol)
    p1r = plot(tSim,rSol[:,1], title="x", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p2r = plot(tSim,rSol[:,2], title="y", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p3r = plot(tSim,rSol[:,3], title="z", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    pr = plot(p1r,p2r,p3r,layout=(3,1))
end

function plotVel(tSim,vSol)
    p1v = plot(tSim,vSol[:,1], title="u", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p2v = plot(tSim,vSol[:,2], title="v", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p3v = plot(tSim,vSol[:,3], title="w", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    pv = plot(p1v,p2v,p3v,layout=(3,1))
end

function plotQuat(tSim,βSol)
    p1 = plot(tSim,βSol[:,1], title="q0", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p2 = plot(tSim,βSol[:,2], title="q1", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p3 = plot(tSim,βSol[:,3], title="q2", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p4 = plot(tSim,βSol[:,4], title="q3", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p = plot(p1,p2,p3,p4,layout=(4,1))
end

function plotQuatDot(tSim,βdotSol)
    p1 = plot(tSim,βdotSol[:,1], title="q0dot", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p2 = plot(tSim,βdotSol[:,2], title="q1dot", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p3 = plot(tSim,βdotSol[:,3], title="q2dot", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p4 = plot(tSim,βdotSol[:,4], title="q3dot", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p = plot(p1,p2,p3,p4,layout=(4,1))
end

function plotAngVel(tSim,ωSol)
    p1 = plot(tSim,ωSol[:,1], title="w1", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p2 = plot(tSim,ωSol[:,2], title="w2", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p3 = plot(tSim,ωSol[:,3], title="w3", xlabel="t")#, label = ["WithUpdate" "WithoutUpdate"], show=true)
    p = plot(p1,p2,p3,layout=(3,1))
end

function plotErrNorm(tSim,βSol)
    errNormβ = zeros(length(tSim),1)
    for i=1:length(tSim)
        errNormβ[i] = norm(βSol[i,:],2)^2 - 1.0
    end
    errNorm = broadcast(abs,errNormβ)
    p = plot(tSim,errNorm, title = "Drift in ||q||")
end
