function adaptGeometry!(p::femProblem,vid::Array{Int64,1},c::Array{Float64,2})
    p.mesh.geometry.coordinates[:,vid]=c;
    b=getBoundary(p.mesh);

    for k in keys(p.degF)
            p.degF[k]=degF(p.mesh, k, b, p.boundary,p.kubPoints);
    end

    if !isempty(p.massM)
            assembMass!(p);
    end
    if !isempty(p.stiffM)
            assembStiff!(p);
    end

    return nothing;
end

function adaptGeometry!(p::femProblem,r::Float64,adaptBoundary::Bool)
        coord=p.mesh.geometry.coordinates;
        dim=p.mesh.topology.D;
        off=p.mesh.topology.offset["$(dim-1)$dim"];
        inc=p.mesh.topology.incidence["$(dim-1)$dim"];

        if adaptBoundary==false
            offe=p.mesh.topology.offset["$(dim-1)0"];
            ince=p.mesh.topology.incidence["$(dim-1)0"];

            bv=Set{Int64}()
            b=Dict{Int64, Array{Int64,1}}();
            for i in 1:(length(off)-1)
                if off[i+1]-off[i]==1
                    incv=ince[offe[i]:(offe[i+1]-1)];
                    for j in incv
                        push!(bv, j);
                    end
                    fi=inc[off[i]];
                    if haskey(b, fi)
                        push!(b[fi],i);
                    else
                        b[fi]=[i];
                    end
                end
            end
            for i in 1:size(coord,2)
                if !in(i,bv)
                    x=r*rand([-1,1])*rand(1)[1];
                    y=rand([-1,1])*sqrt(r*r-x*x)*rand(1)[1];
                    coord[:,i]=[coord[1,i]+x, coord[2,i]+y];
                end
            end
        else
            b=getBoundary(p.mesh);
            for i in 1:size(coord,2)
                x=r*rand([-1,1])*rand(1)[1];
                y=rand([-1,1])*sqrt(r*r-x*x)*rand(1)[1];
                coord[:,i]=[coord[1,i]+x, coord[2,i]+y];
            end
        end

        for k in keys(p.degF)
            p.degF[k]=degF(p.mesh, k, b, p.boundary,p.kubPoints);
        end
        if !isempty(p.massM)
                assembMass!(p);
        end
        if !isempty(p.stiffM)
                assembStiff!(p);
        end
        return nothing;
end

function adaptGeometry!(p::femProblem,pertX::Float64, pertY::Float64,adaptBoundary::Bool)
        coord=p.mesh.geometry.coordinates;
        dim=p.mesh.topology.D;
        off=p.mesh.topology.offset["$(dim-1)$dim"];
        inc=p.mesh.topology.incidence["$(dim-1)$dim"];
        xR=p.mesh.geometry.r[1]; xL=p.mesh.geometry.l[1]; yR=p.mesh.geometry.r[2]; yL=p.mesh.geometry.l[2]
        dx=(xR-xL)/p.mesh.topology.n[1];
        dy=(yR-yL)/p.mesh.topology.n[2];

        if adaptBoundary==false
            offe=p.mesh.topology.offset["$(dim-1)0"];
            ince=p.mesh.topology.incidence["$(dim-1)0"];

            bv=Set{Int64}()
            b=Dict{Int64, Array{Int64,1}}();
            for i in 1:(length(off)-1)
                if off[i+1]-off[i]==1
                    incv=ince[offe[i]:(offe[i+1]-1)];
                    for j in incv
                        push!(bv, j);
                    end
                    fi=inc[off[i]];
                    if haskey(b, fi)
                        push!(b[fi],i);
                    else
                        b[fi]=[i];
                    end
                end
            end
            for i in 1:size(coord,2)
                if !in(i,bv)
                    x=(coord[1,i]);
                    y=(coord[2,i]);
                    hx=coord[1,i]+pertX*sin(2*pi*(y-yL)/(yR-yL))*dx;
                    hy=coord[2,i]+pertY*sin(2*pi*(x-xL)/(xR-xL))*dy;
                    coord[:,i]=[hx, hy];
                end
            end
        else
            b=getBoundary(p.mesh);
            for i in 1:size(coord,2)
                x=(coord[1,i]);
                y=(coord[2,i]);
                hx=coord[1,i]+pertX*sin(2*pi*(y-yL)/(yR-yL))*dx;
                hy=coord[2,i]+pertY*sin(2*pi*(x-xL)/(xR-xL))*dy;
                coord[:,i]=[hx, hy];
            end
        end

        for k in keys(p.degF)
            p.degF[k]=degF(p.mesh, k, b, p.boundary,p.kubPoints);
        end
        if !isempty(p.massM)
                assembMass!(p);
        end
        if !isempty(p.stiffM)
                assembStiff!(p);
        end
        return nothing;
end
