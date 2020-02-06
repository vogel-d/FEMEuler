function plotSolution(p::femProblem, key::Symbol, t::Float64, comp::Int=0;filename::String="", c::Symbol=:plasma)
    nf=p.mesh.topology.size[3];

    off=p.mesh.topology.offset["20"];
    inc=p.mesh.topology.incidence["20"];

    sol=getfield(p.solution[t],key);
    skm=getSolution(sol,p.degFBoundary[p.femType[key][1]],p.mesh,p.femType[key][1],comp)

    pc=plot();
    for i in 0:0.001:1
        color=ColorGradient(c)[i];
        plot!([i,i],[0,1],linecolor=color, legend=false, yaxis=nothing, xaxis=nothing, foreground_color=:white)
    end

    smin=minimum(skm);
    smax=maximum(skm);
    mw=0.5*(smin+smax);
    mr=0.5*(mw+smax);
    ml=0.5*(smin+mw);
    a=[(0.0,1.0,text("$(round(smin,digits=3))",:bottom, 6)),(0.25,1.0,text("$(round(ml,digits=3))",:bottom,6)),(0.5,1.0,text("$(round(mw,digits=3))",:bottom,6)),(0.75,1.0,text("$(round(mr,digits=3))",:bottom,6)),(1.0,1.0,text("$(round(smax,digits=3))",:bottom,6))];
    annotate!(a)

    pl=plot();
    for k in 1:nf
        fillPoly(p.mesh.geometry.coordinates[:,inc[off[k]:off[k+1]-1]], ColorGradient(c)[(skm[k]-smin)/(smax-smin)]);
    end

    pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
    if filename==""
        display(pf)
    else
        savefig(pf, filename);
    end
end

function getSolution(sol::Array{Float64,1},degF::degF{1}, m::mesh, key::Symbol, comp::Int)
    nf=m.topology.size[3];
    skm=zeros(Float64,nf);
    fComp=getElementProperties(key,m.meshType,0.5,0.5);
    for k in 1:nf
        sk=sol[l2g(degF, k)];
        skm[k]=dot(fComp,sk);
    end
    return skm;
end

function getSolution(sol::Array{Float64,1},degF::degF{2}, m::mesh, key::Symbol, comp::Int)
    comp==0 && error("Für eine vektorielle Variable muss die zu plottende Komponente spezifiziert werden, für x comp=1 und für y comp=2.")
    nf=m.topology.size[3];
    skm=zeros(Float64,nf);
    J=Array{Float64,2}(undef,2,2);
    dJ=0.0;
    coord=Array{Float64,2}(undef,2,m.meshType);
    fComp=getElementProperties(key,m.meshType,0.5,0.5);
    for k in 1:nf
        sk=sol[l2g(degF, k)];
        dJ=jacobi!(J,m,k,0.5,0.5,coord);
        fLoc=((1/dJ)*J*fComp)*sk
        skm[k]=fLoc[comp];
    end
    return skm;
end

function fillPoly(c::Array{Float64,2}, fc::RGBA{Float64}; fillrange::Int64=0, fillalpha::Float64=1.0)
    x=push!(c[1,:],c[1,1]);
    y=push!(c[2,:],c[2,1]);
    p=Plots.plot!(x,y,fill=(fillrange,fillalpha,fc),legend=false, linecolor=fc, linewidth=0.0);
    return(p)
end
