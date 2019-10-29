function generatePeriodicBoundary!(p::femProblem)
    dF=deepcopy(p.degF);
    equals=p.equals;
    #periodische Randbedingung:
    #über equals wird erfasst, welches Randfacet mit welchem Randfacet "gleichzusetzen" ist,
    #dann wird über die Inzidenz Kante->Freiheitsgrad erfasst welche Freiheitsgrade mit welchen Freiheitsgraden
    #gleichzusetzen sind und die Inzidenz angepasst.
    for k in keys(p.degF)
        lc=dF[k].num;
        coord=dF[k].coordinates[:,1:lc];
        bcoord=dF[k].coordinates[:,(lc+1):end];
        lmax=size(dF[k].coordinates,2);
        inc=dF[k].incidence;
        incn=spzeros(length(inc))
        linc=length(inc);
        isempty(bcoord) && continue;
        done=Set{Int64}();
        eb=dF[k].edgeBoundaryIncidence;
        be=dF[k].boundaryEdgeIncidence;
        for j in 1:length(be)
            edges=be[j];
            for e1 in edges
                e2=equals[e1];
                e2==0 && continue;
                b1=eb[e1];
                b2=eb[e2];

                hp=bcoord[:,b1].==bcoord[:,b2]
                while sum(hp)==0
                    pushfirst!(b2,pop!(b2));
                    hp=bcoord[:,b1].==bcoord[:,b2]
                end

                for i in 1:length(b1)
                    if !in(b1[i],done) && !in(b2[i],done)
                        coord=hcat(coord,bcoord[:,b1[i]]);

                        f=findall(b1[i]+lc, inc);
                        append!(f,findall(b2[i]+lc, inc));
                        incn[f].=size(coord,2);

                        push!(done, b1[i]);
                        push!(done, b2[i]);
                    end
                end
            end
        end
        for i in 1:length(incn)
            if incn[i]==0
                incn[i]=inc[i];
            end
        end
        bcoord = bcoord[:,setdiff(1:size(bcoord,2), done)];

        xR=p.mesh.geometry.r[1]; xL=p.mesh.geometry.l[1]; yR=p.mesh.geometry.r[2]; yL=p.mesh.geometry.l[2]

        pos=findall([xL,yL],coord);
        missing=Set{Int64}();
        for i in 1:size(coord,2)
            if isapprox(coord[:,i],[xL,yR])
                push!(missing,i);
                replace!(incn, i=>pos[1]);
            elseif isapprox(coord[:,i],[xR,yL])
                push!(missing,i);
                coord[:,i]=[xL,yL];
                replace!(incn, i=>pos[1]);
            elseif isapprox(coord[:,i],[xR,yR])
                push!(missing,i);
                coord[:,i]=[xL,yL];
                replace!(incn, i=>pos[1]);
            end
        end

        range=Set(1:size(coord,2));
        appeared=sort(collect(setdiff(range,missing)));

        coord=coord[:,sort(collect(appeared))];
        for j in 1:length(incn)
            for i in missing
                if incn[j]>i
                    incn[j]-=1;
                end
            end
            if incn[j]>size(coord,2)
                incn[j]=pos[1];
            end
        end

        dF[k].coordinates=coord;
        dF[k].incidence=incn;
        dF[k].num=size(coord,2);
    end
    p.degFBoundary=dF;
end
