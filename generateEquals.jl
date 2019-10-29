function generateEquals!(p::femProblem, cond::Tuple{Symbol,Symbol})
    p.boundaryCondition=cond;
    m=p.mesh;
    type=m.meshType;
    nx=m.topology.n[1];
    ny=m.topology.n[2];
    ne=m.topology.size[2];
    equals=spzeros(Int64, ne);
    if type==4
        if cond==(:periodic, :periodic)
            for i in 1:nx
                equals[i]=i+nx*ny;
            end
            for i in 1:ny
                equals[i*(nx+1)+nx*ny]=i*(nx+1)+nx*(ny+1);
            end
        elseif cond==(:periodic, :constant)
            for i in 1:ny
                equals[i*(nx+1)+nx*ny]=i*(nx+1)+nx*(ny+1);
            end
        elseif cond==(:constant, :periodic)
            for i in 1:nx
                equals[i]=i+nx*ny;
            end
        end
    elseif type==3
        if isodd(ny)
            ns==true && error("Für ein Dreiecksmesh mit ungerader Anzahl an Elementen in y-Richtung können die periodischen Ranbedingung
                                nicht auf die oberen und unteren Kanten des Meshes angewendet werden, da die Anzahl der Randfreiheitsgrade nicht übereinstimmt.");
            a=Int64(0.5*(ny+1)*(2*nx-1));
            for i in 1:ny
                equals[a+1+(i-1)*2*nx]=a+i*2*nx;
            end
        else
            ny2=Int64(0.5*ny);
            if cond==(:periodic, :periodic)
                for i in 1:nx
                    equals[i]=i+ny2*(2*nx-1);
                end
                a=(ny2+1)*nx+ny2*(nx-1)
                for i in 1:ny
                    equals[a+1+(i-1)*2*nx]=a+i*2*nx;
                end
            elseif cond==(:periodic, :constant)
                a=(ny2+1)*nx+ny2*(nx-1)
                for i in 1:ny
                    equals[a+1+(i-1)*2*nx]=a+i*2*nx;
                end
            elseif cond==(:constant, :periodic)
                for i in 1:nx
                    equals[i]=i+ny2*(2*nx-1);
                end
            end
        end
    end
    p.equals=equals;
    return nothing;
end
