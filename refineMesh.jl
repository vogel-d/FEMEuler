#rx, /ry ist die Anzahl der Unterteilung einer Zelle des meshes m in x/y-Richtung
function refineRectMesh(m::mesh,rx::Int,ry::Int,condEW::Symbol, condTB::Symbol)
    nx=m.topology.n[1]; ny=m.topology.n[2]
    inc=Int[];
    off=collect(1:rx*ry:m.topology.size[3]*rx*ry+1) #statt 4 n?
    z=1
    for i in 1:ny
        for j in 1:nx
            for k in 0:nx*rx:nx*rx*(ry-1)
                for l in 0:rx-1
                    push!(inc,z+k+l);
                end
            end
            z+=rx;
        end
        z+=rx*(ry-1)*nx;
    end
    mf=generateRectMesh(rx*nx,ry*ny,condEW, condTB, m.geometry.l[1], m.geometry.r[1], m.geometry.l[2], m.geometry.r[2])
    mf.topology.incidence["CF"]=inc #CF = Coarse -> Fine
    mf.topology.offset["CF"]=off
    return mf;
end

function refineCubedSphere(m::mesh,rn::Int)
    n=m.topology.n[1];
    inc=Int[];
    off=collect(1:rn^2:m.topology.size[3]*rn^2+1) #statt 4 n?
    z=1
    for h in 1:6
        for i in 1:n
            for j in 1:n
                for k in 0:n*rn:n*rn*(rn-1)
                    for l in 0:rn-1
                        push!(inc,z+k+l);
                    end
                end
                z+=rn;
            end
            z+=rn*(rn-1)*n;
        end
    end
    mf=generateCubedSphere(rn*n,m.geometry.r[1])
    mf.topology.incidence["CF"]=inc #CF = Coarse -> Fine
    mf.topology.offset["CF"]=off
    return mf;
end
