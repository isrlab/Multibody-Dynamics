using PyPlot;

function plotPos(tSim,rSol, titleStr::String="Pos")
    figure(); clf;
    pygui(true);
    subplot(3,1,1);
    PyPlot.plot(tSim,rSol[:,1]); xlabel("t"); ylabel("x");
    subplot(3,1,2);
    PyPlot.plot(tSim,rSol[:,2]); xlabel("t"); ylabel("y")
    subplot(3,1,3);
    PyPlot.plot(tSim,rSol[:,3]); xlabel("t"); ylabel("z")
    PyPlot.suptitle(titleStr)
    tight_layout();
end

function plotVel(tSim,vSol, titleStr::String="Vel")
    figure(); clf;
    pygui(true);
    subplot(3,1,1);
    PyPlot.plot(tSim,vSol[:,1]); xlabel("t"); ylabel("u");
    subplot(3,1,2);
    PyPlot.plot(tSim,vSol[:,2]); xlabel("t"); ylabel("v")
    subplot(3,1,3);
    PyPlot.plot(tSim,vSol[:,3]); xlabel("t"); ylabel("w")
    PyPlot.suptitle(titleStr)
    tight_layout();

end

function plotQuat(tSim,βSol, titleStr::String="Quaternion")
    figure(); clf;
    pygui(true);
    subplot(4,1,1);
    PyPlot.plot(tSim,βSol[:,1]); xlabel("t"); ylabel("β0");
    subplot(4,1,2);
    PyPlot.plot(tSim,βSol[:,2]); xlabel("t"); ylabel("β1")
    subplot(4,1,3);
    PyPlot.plot(tSim,βSol[:,3]); xlabel("t"); ylabel("β2")
    subplot(4,1,4);
    PyPlot.plot(tSim,βSol[:,4]); xlabel("t"); ylabel("β3")
    PyPlot.suptitle(titleStr)
    tight_layout();

end

function plotQuatDot(tSim,βdotSol, titleStr::String="Quat Dot")
    figure(); clf;
    pygui(true);
    subplot(4,1,1);
    PyPlot.plot(tSim,βdotSol[:,1]); xlabel("t"); ylabel("̇β̇0");
    subplot(4,1,2);
    PyPlot.plot(tSim,βdotSol[:,2]); xlabel("t"); ylabel("β̇1")
    subplot(4,1,3);
    PyPlot.plot(tSim,βdotSol[:,3]); xlabel("t"); ylabel("β̇2")
    subplot(4,1,4);
    PyPlot.plot(tSim,βdotSol[:,4]); xlabel("t"); ylabel("β̇3")
    PyPlot.suptitle(titleStr)
    tight_layout();

end

function plotAngVel(tSim,ωSol, titleStr::String="Ang Vel")
    figure(); clf;
    pygui(true);
    subplot(3,1,1);
    PyPlot.plot(tSim,ωSol[:,1]); xlabel("t"); ylabel("ω1");
    subplot(3,1,2);
    PyPlot.plot(tSim,ωSol[:,2]); xlabel("t"); ylabel("ω2")
    subplot(3,1,3);
    PyPlot.plot(tSim,ωSol[:,3]); xlabel("t"); ylabel("ω3")
    PyPlot.suptitle(titleStr)
    tight_layout();

end

function plotErrNorm(tSim,βSol)
    errNormβ = zeros(length(tSim),1)
    for i=1:length(tSim)
        errNormβ[i] = norm(βSol[i,:],2)^2 - 1.0
    end
    errNorm = broadcast(abs,errNormβ)
    figure(); clf;
    pygui(true)
    plot(tSim,errNorm); title("Drift in ||q||")
    tight_layout();

end
