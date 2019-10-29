function generateMixedBoundary!(p::femProblem)
    dF=deepcopy(p.degF);
    equals=p.equals;
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
        bcoord = bcoord[:,setdiff(1:size(bcoord,2), done)];
        bcoord2=p.degF[k].coordinates[:,(p.degF[k].num+1):end]
        for i in 1:length(incn)
            if incn[i]==0
                if inc[i]<=lc
                    incn[i]=inc[i];
                else
                    b=findall(bcoord2[:,inc[i]-lc],bcoord);
                    incn[i]=b[1]+size(coord,2);
                end
            end
        end

        dF[k].coordinates=hcat(coord,bcoord);
        dF[k].incidence=incn;
        dF[k].num=size(coord,2);
    end
    p.degFBoundary=dF;
end
