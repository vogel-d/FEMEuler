mutable struct degF{N}
    coordinates::Array{Float64,2};
    num::Int64;
    edgeBoundaryIncidence::SparseVector{Array{Int64,1},Int64};
    boundaryEdgeIncidence::Array{Array{Int64,1},1};
    incidence::Array{Int64,1};
    offset::Array{Int64,1};
    referenceCoordinates::Array{Float64,2};
    referenceBoundary::Array{Float64,2};
    phi::Array{Float64,N};
    divphi::Array{Float64,3};
    gradphi::Array{Float64,4};
    discont::Bool;
    components::Array{Int64,1};
end


function degF(m::mesh, femType::Symbol, boundary::Dict{Int64, Array{Int64,1}}, fset::Set{Int64}, kubPoints::Array{Float64,2})
    D=m.topology.D;
    a=0;
    nf=m.topology.size[D+1];

    if m.topology.offset["$D$a"][2]-m.topology.offset["$D$a"][1]==4
        phi, divphi,  gradphi, cref,ccenters, comp, discontType=getQuadElementProperties(femType, kubPoints)
    else
        phi, divphi,  gradphi, cref,ccenters, comp, discontType=getTriElementProperties(femType, kubPoints)
    end

    degComp=Int64[];
    degCompB=Int64[];
    trans(co,x,y)=transformation(m,co,x,y);
    degFc=Array{Float64,2}(undef,2,0);
    degFbc=Array{Float64,2}(undef,2,0);
    degFi=Int64[];
    degFo=Array{Int64}(1:size(cref,2):(size(cref,2)*nf+1));
    degFbe=Array{Array{Int64,1},1}();
    degFeb=spzeros(Array{Int64,1}, m.topology.size[2]);
    h=Int64[];
    z1=false;
    z2=false;

    if discontType
        z3=false;
        z4=false;
    else
        z3=true;
        z4=true;
    end
    #Konzept: Herausfinden der am Rand liegenden Freiheitsgrade
    #durch Zuordnung der Freiheitsgrade an Kanten im Referenzelement
    #(Die Randkanten liefert boundary)
    for k in 1:nf
        rstart=m.topology.offset["$D$a"][k];
        rend=m.topology.offset["$D$a"][k+1]-1;
        vertices=m.topology.incidence["$D$a"][rstart:rend];
        coord=m.geometry.coordinates[:,vertices];
        if in(k,fset) #falls das Zelle ein Randfacet ist
            incz=boundary[k]; #werden Randkanten der Zelle bestimmt
            bdegF=Int64[];
            bdegFh=Int64[];
            offb=Int64[1];
            zb=1;
            for i in 1:length(incz) #und dann für jede Randkante
                rstart=m.topology.offset["10"][incz[i]];
                cz=m.geometry.coordinates[:,m.topology.incidence["10"][rstart:(rstart+1)]];
                mz=(1/size(cz,2)).*[sum(cz[1,:]), sum(cz[2,:])]; #der Mittelpunkt bestimmt
                for j in 1:size(ccenters,1)
                    #danach wird geprüft welcher Referenzmittelpunkt nach der Transformation
                    #diesen Mittelpunkt ergibt
                    if isapprox(mz,trans(coord,ccenters[j,1],ccenters[j,2]))
                        for l in 3:size(ccenters,2)
                            if ccenters[j,l]==1.0
                                push!(bdegF,l-2);
                                push!(bdegFh,incz[i]);
                                zb+=1;
                                #dann werden alle Freiheitsgrade die diesem
                                #Referenzmittelpunkt zugeordnet sind erfasst
                            end
                        end
                    end
                end
                push!(offb,zb);
            end
            #dann werden alle untersucht
            bdegFg=spzeros(Int64, length(bdegF));
            for i in 1:size(cref,2)
                if in(i,bdegF)
                    #falls es ein Randfreiheitsgrad ist wird er mit höheren Werten markiert
                    t=trans(coord,cref[1,i],cref[2,i]);
                    tcomp=comp[i];
                    if z1
                        b=findall(t,degFbc);
                        bh=true;
                        if length(b)>0
                            for ic in 1:length(b)
                                if tcomp==degCompB[b[ic]]
                                    push!(degFi,b[ic]);
                                    push!(h, length(degFi));
                                    ind=findall(i,bdegF);
                                    append!(degFbe[b[ic]],bdegFh[ind]);
                                    bdegFg[ind].=b[ic];
                                    bh=false;
                                    break;
                                end
                            end
                        end
                        if length(b)==0 || bh
                            degFbc=hcat(degFbc,[t[1], t[2]]);
                            push!(degFi, size(degFbc,2));
                            push!(degCompB, tcomp);
                            push!(h, length(degFi));
                            ind=findall(i,bdegF);
                            push!(degFbe,bdegFh[ind]);
                            bdegFg[ind].=size(degFbc,2);
                        end
                    else
                        degFbc=hcat(degFbc,[t[1], t[2]]);
                        push!(degFi, size(degFbc,2));
                        push!(degCompB, tcomp);
                        push!(h, length(degFi));
                        ind=findall(i,bdegF);
                        push!(degFbe,bdegFh[ind]);
                        bdegFg[ind].=size(degFbc,2);
                        z1=z3;
                    end
                else
                    #falls nicht findet die normale Erzeugung als Freiheitsgrad statt
                    t=trans(coord,cref[1,i],cref[2,i]);
                    tcomp=comp[i];
                    if z2
                        b=findall(t,degFc);
                        bh=true;
                        if length(b)>0
                            for ic in 1:length(b)
                                if tcomp==degComp[b[ic]]
                                    push!(degFi,b[ic])
                                    bh=false;
                                    break;
                                end
                            end
                        end
                        if length(b)==0 || bh
                            degFc=hcat(degFc, [t[1], t[2]]);
                            push!(degFi, size(degFc,2));
                            push!(degComp,tcomp)
                        end
                    else
                        degFc=hcat(degFc, [t[1], t[2]]);
                        push!(degFi, size(degFc,2));
                        push!(degComp,tcomp)
                        z2=z4;
                    end
                end
            end
            for i in 1: length(incz)
                degFeb[incz[i]]=bdegFg[offb[i]:(offb[i+1]-1)];
            end
        else
            #falls das Facet nicht zum Rand gehört
            #werden alle Freiheitsgrade normal erzeugt
            for i in 1:size(cref,2)
                t=trans(coord,cref[1,i],cref[2,i]);
                tcomp=comp[i];
                if z2
                    b=findall(t,degFc);
                    bh=true;
                    if length(b)>0
                        for ic in 1:length(b)
                            if tcomp==degComp[b[ic]]
                                push!(degFi,b[ic])
                                bh=false;
                                break;
                            end
                        end
                    end
                    if length(b)==0 || bh
                        degFc=hcat(degFc, [t[1], t[2]]);
                        push!(degFi, size(degFc,2));
                        push!(degComp,tcomp)
                    end
                else
                    degFc=hcat(degFc, [t[1], t[2]]);
                    push!(degFi, size(degFc,2));
                    push!(degComp,tcomp)
                end
            end
        end
    end
    s=size(degFc,2);
    for i in 1:length(h)
        degFi[h[i]]+=s;
    end
    lc=size(degFc,2);
    degFc=hcat(degFc,degFbc)
    degF(degFc, lc, degFeb, degFbe, degFi, degFo, cref, ccenters, phi, divphi, gradphi, discontType, comp);
end
