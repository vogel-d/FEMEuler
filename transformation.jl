function transformation!(r::Array{Float64,1}, m::mesh, coord, x::Float64, y::Float64)
    #es wird von einem Gitter mit einheitlichen Elementen ausgegangen
    #davon ausgehend lässt sich n allgemein durch die ersten beiden Offset-Einträge feststellen
    n=m.meshType;
    if n==3
        for i in 1:length(r)
            r[i]=coord[i,1]+(coord[i,2]-coord[i,1])*x+(coord[i,3]-coord[i,1])*y
        end
        return nothing
    elseif n==4
        if m.geometry.dim==2
            for i in 1:length(r)
                r[i]=coord[i,1]+(coord[i,2]-coord[i,1])*x+(coord[i,4]-coord[i,1])*y+(coord[i,3]-coord[i,4]-coord[i,2]+coord[i,1])*x*y;
            end
            return nothing
        else
            k=copy(r)
            for i in 1:length(r)
                r[i]=coord[i,1]+(coord[i,2]-coord[i,1])*x+(coord[i,4]-coord[i,1])*y+(coord[i,3]-coord[i,4]-coord[i,2]+coord[i,1])*x*y;
            end
            r[:]=m.geometry.r[1].*r./norm(r,2);
            return nothing
        end
    else
        error("Diese Funktion benötigt ein Mesh mit drei- oder viereckigen Elementen");
    end
end

function transformation(m::mesh, coord, x::Float64, y::Float64)
    #es wird von einem Gitter mit einheitlichen Elementen ausgegangen
    #davon ausgehend lässt sich n allgemein durch die ersten beiden Offset-Einträge feststellen
    n=m.meshType;
    if n==3
        r1=coord[:,1];
        return r1+(coord[:,2]-r1)*x+(coord[:,3]-r1)*y
    elseif n==4
        r1=coord[:,1];
        if m.geometry.dim==2
            return r1+(coord[:,2]-r1)*x+(coord[:,4]-r1)*y+(coord[:,3]-coord[:,4]-coord[:,2]+r1)*x*y
        else
            r=r1+(coord[:,2]-r1)*x+(coord[:,4]-r1)*y+(coord[:,3]-coord[:,4]-coord[:,2]+r1)*x*y;
            return m.geometry.r[1]*r/norm(r,2);
        end
    else
        error("Diese Funktion benötigt ein Mesh mit drei- oder viereckigen Elementen");
    end
end
