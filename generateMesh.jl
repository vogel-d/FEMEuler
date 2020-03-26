#Funktion zum Generieren von zweidimensionalen Rechteck-Gittern
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       xl bzw. yl ist die Länge des Meshes in x- bzw. y-Richtung, als Default ist die Länge nx bzw. ny
function generateRectMesh(nx::Int, ny::Int, condEW::Symbol, condTB::Symbol, xl::Float64=0.0, xr::Float64=Float64(nx), yl::Float64=0.0, yr::Float64=Float64(ny))

    #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
    size=[(ny+1)*(nx+1), ny*(nx+1)+nx*(ny+1), nx*ny];
    nk=4;

    #Initialisieren des Offsets mit den Einträgen "20" und "10"
    offe=collect(1:2:(2*size[2]+1));
    off=Dict("20"=>collect(1:nk:(nk*size[3]+1)),"10"=>offe);

    #Berechnen der Koordinatenmatrix, basierend auf einem äquidistanten Gitter
    coord=Array{Float64}(undef,2,size[1]);
    ax=(xr-xl)/nx;
    ay=(yr-yl)/ny;
    zx=xl;
    zy=yl;
    z=1;
    for k in 1:(ny+1)
        for l in 1:(nx+1)
            coord[:,z]=[zx, zy];
            zx+=ax;
            z+=1;
        end
        zx=xl;
        zy+=ay;
    end

    #Berechnen der Inzidenz 2->0
    incf=Int[];
    z=1;
    for k in 1:ny
        for h in 1:nx
            i=[z, z+1, z+1+nx+1, z+nx+1];
            append!(incf,i);
            z+=1;
        end
        z+=1;
    end

    #Berechnen der Inzidenz 1->0
    ince=Int[];
    z=1;
    for k in 1:(nx+1)*ny
        i=[z, z+nx+1];
        append!(ince,i);
        z+=1;
    end
    zh=1;
    for h in 1:nx
        z=zh;
        for k in 1:(ny+1)
            i=[z, z+1]
            append!(ince,i);
            z+=nx+1
        end
        zh+=1;
    end

    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>incf,"10"=>ince);

    #Initialisieren des Boundaryvektoren mit Einträgen
    #   0 für innere Kante/Knoten,
    #   1 für RandKante/Knoten mit freeslip,
    #   -x für periodische RandKante/Knoten mit x gegenüberliegende RandKante/Knoten
    bE=spzeros(Int, size[2]);
    bV=spzeros(Int, size[1]);
    ordV=Int[];
    corners=Set{Int}();
    #Berechnen des Boundaryvektors für die Randknoten
    for i in 1:size[1]
        if (coord[1,i]==xl || coord[1,i]==xr) && (coord[2,i]==yl || coord[2,i]==yr)
            push!(corners,i);
            if condEW==:constant && condTB==:constant
                bV[i]=1;
            elseif condEW==:periodic && condTB==:periodic && !(coord[1,i]==xl && coord[2,i]==yl)
                bV[i]=-1;
            elseif condEW==:periodic && condTB==:constant && coord[1,i]==xl
                bV[i]=-i-nx;
            elseif condTB==:periodic && condEW==:constant && coord[2,i]==yl
                bV[i]=-i-(nx+1)*ny;
            end
        elseif condEW==:constant && (coord[1,i]==xl || coord[1,i]==xr)
            bV[i]=1;
        elseif coord[1,i]==xl
            #condEW==:periodic
            bV[i]=-i-nx;
        elseif condTB==:constant && (coord[2,i]==yl || coord[2,i]==yr)
            bV[i]=1;
        elseif coord[2,i]==yl
            #condTB==:periodic
            bV[i]=-i-(nx+1)*ny;
        end
    end

    #Berechnen des Boundaryvektors für die Randkanten
    for i in 1:size[2]
        v1=ince[offe[i]];
        v2=ince[offe[i]+1]
        b1=bV[v1];
        b2=bV[v2];
        if b1==0 && b2==0
            continue;
        elseif v1==1 || v2==1 #Kante unten links gesondert behandelt wegen periodic/periodic, 0 in bV[1] macht Probleme
            if b2==b1==1
                bE[i]=1;
            elseif (b1==1 && b2<=0) || (b1<=0 && b2==1)
                bE[i]=1;
            elseif (b1==0 && b2<0) || (b1<0 && b2==0) || (b1<0 && b2<0)
                if coord[1,v1]==coord[1,v2]
                    bE[i]=-i-nx;
                elseif coord[2,v1]==coord[2,v2]
                    bE[i]=-i-ny;
                end
            end
        elseif (b1!=0 && b2!=0) || in(v1,corners) || in(v2,corners)
            if b2==b1==1
                bE[i]=1;
            elseif (b1==1 && b2<=0) || (b1<=0 && b2==1)
                bE[i]=1;
            elseif b1<0 && b2<0
                if coord[1,v1]==coord[1,v2]
                    bE[i]=-i-nx;
                elseif coord[2,v1]==coord[2,v2]
                    bE[i]=-i-ny;
                end
            end
        end
    end


    #Initialisieren der Topologie, Geometrie und damit des Meshes
    n=Int[nx,ny];
    l=Float64[xl,yl];
    r=Float64[xr,yr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV);

    return m
end

#Funktion zum Generieren von zweidimensionalen Dreieck-Gitter im Rahmen eines Rechtecks
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       xl bzw. yl ist die Länge des Meshes in x- bzw. y-Richtung, als Default ist die Länge nx bzw. ny
function generateTriMesh(nx::Int, ny::Int, condEW::Symbol, condTB::Symbol, xl::Float64=0.0, xr::Float64=Float64(nx), yl::Float64=0.0, yr::Float64=Float64(ny))

    #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
    size=[(ny+1)*(nx+1), ny*(nx+1)+nx*(ny+1)+nx*ny, 2*nx*ny];
    nk=3;

    #Initialisieren des Offsets mit den Einträgen "20" und "10"
    offe=collect(1:2:(2*size[2]+1));
    off=Dict("20"=>collect(1:nk:(nk*size[3]+1)),"10"=>offe);

    #Berechnen der Koordinatenmatrix, basierend auf einem äquidistanten Gitter
    coord=Array{Float64}(undef,2,size[1]);
    ax=(xr-xl)/nx;
    ay=(yr-yl)/ny;
    z=1;
    for k in 1:(ny+1)
        for l in 1:(nx+1)
            coord[:,z]=[xl+(l-1)*ax,yl+(k-1)*ay]
            z+=1;
        end
    end

    #Berechnen der Inzidenz 2->0
    incf=Int[];
    z=1;
    for k in 1:ny
        for h in 1:nx
            i=[z, z+1, z+1+nx];
            append!(incf,i);
            i=[z+1, z+nx+1, z+nx+2]
            append!(incf,i);
            z+=1;
        end
        z+=1;
    end


    #Berechnen der Inzidenz 1->0
    ince=Int[];
    z=1;
    for k in 1:(nx+1)*ny
        i=[z, z+nx+1];
        append!(ince,i);
        z+=1;
    end
    zh=1;
    for h in 1:nx
        z=zh;
        for k in 1:(ny+1)
            i=[z, z+1]
            append!(ince,i);
            z+=nx+1
        end
        zh+=1;
    end

    #schräge Kanten
    for k in 1:(nx+1)*ny
        if mod(k,nx+1)!=0
            i=[k+1, k+nx+1];
            append!(ince,i);
        end
    end
    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>incf,"10"=>ince);


    corners=Set{Int}();
    bE=spzeros(Int, size[2]);
    bV=spzeros(Int, size[1]);
    #berechnen des boundary Vektors für die Randknoten
    for i in 1:size[1]
        if (coord[1,i]==xl || coord[1,i]==xr) && (coord[2,i]==yl || coord[2,i]==yr)
            push!(corners,i);
            if condEW==:constant && condTB==:constant
                bV[i]=1;
            elseif condEW==:periodic && condTB==:periodic && !(coord[1,i]==xl && coord[2,i]==yl)
                bV[i]=-1;
            elseif condEW==:periodic && condTB==:constant && coord[1,i]==xl
                bV[i]=-i-nx;
            elseif condTB==:periodic && condEW==:constant && coord[2,i]==yl
                bV[i]=-i-(nx+1)*ny;
            end
        elseif condEW==:constant && (coord[1,i]==xl || coord[1,i]==xr)
            bV[i]=1;
        elseif coord[1,i]==xl
            #condEW==:periodic
            bV[i]=-i-nx;
        elseif condTB==:constant && (coord[2,i]==yl || coord[2,i]==yr)
            bV[i]=1;
        elseif coord[2,i]==yl
            #condTB==:periodic
            bV[i]=-i-(nx+1)*ny;
        end
    end

    #berechnen des boundary Vektors für die Randkanten
    for i in 1:size[2]
        v1=ince[offe[i]];
        v2=ince[offe[i]+1]
        b1=bV[v1];
        b2=bV[v2];
        if b1==0 && b2==0
            continue;
        elseif v1==1 || v2==1 #Kante unten links gesondert behandelt wegen periodic/periodic, 0 in bV[1] macht Probleme
            if b2==b1==1
                bE[i]=1;
            elseif (b1==1 && b2<=0) || (b1<=0 && b2==1)
                bE[i]=1;
            elseif (b1==0 && b2<0) || (b1<0 && b2==0) || (b1<0 && b2<0)
                if coord[1,v1]==coord[1,v2]
                    bE[i]=-i-nx;
                elseif coord[2,v1]==coord[2,v2]
                    bE[i]=-i-ny;
                end
            end
        elseif (b1!=0 && b2!=0) || in(v1,corners) || in(v2,corners)
            if b2==b1==1
                bE[i]=1;
            elseif (b1==1 && b2<=0) || (b1<=0 && b2==1)
                if !(coord[1,v1]!=coord[1,v2] && coord[2,v1]!=coord[2,v2]) #nur Kante die tatsächlich am Rand liegt (keine schrägen)
                    bE[i]=1;
                end
            elseif b1<0 && b2<0
                if coord[1,v1]==coord[1,v2]
                    bE[i]=-i-nx;
                elseif coord[2,v1]==coord[2,v2]
                    bE[i]=-i-ny;
                end
            end
        end
    end

    #Initialisieren der Topologie, Geometrie und damit des Meshes
    n=Int[nx,ny];
    l=Float64[xl,yl];
    r=Float64[xr,yr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV);

    return m
end


function generateTriMesh2(nx::Int, ny::Int, condEW::Symbol, condTB::Symbol, xl::Float64=0.0, xr::Float64=Float64(nx), yl::Float64=0.0, yr::Float64=Float64(ny))

    #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
    size=[(ny+1)*(nx+1)+nx*ny, ny*(nx+1)+nx*(ny+1)+4*nx*ny, 4*nx*ny];
    nk=3;

    #Initialisieren des Offsets mit den Einträgen "20" und "10"
    offe=collect(1:2:(2*size[2]+1));
    off=Dict("20"=>collect(1:nk:(nk*size[3]+1)),"10"=>offe);

    #Berechnen der Koordinatenmatrix, basierend auf einem äquidistanten Gitter
    coord=Array{Float64}(undef,2,size[1]);
    ax=(xr-xl)/nx;
    ay=(yr-yl)/ny;
    z=1;
    for k in 1:(ny+1)
        for l in 1:(nx+1)
            coord[:,z]=[xl+(l-1)*ax,yl+(k-1)*ay]
            z+=1;
        end
    end
    for k in 1:ny
        for l in 1:nx
            coord[:,z]=[xl+(l-1+0.5)*ax,yl+(k-1+0.5)*ay]
            z+=1
        end
    end

    #Berechnen der Inzidenz 2->0
    incf=Int[];
    z=1;
    for k in 1:ny
        for l in 1:nx
            #viereck in vier dreiecke zerteilt (über diagonalen)
            #   |\/|    = eine Viereck-Zelle
            #   |/\|    = vier Dreieck-Zellen

            #vertex-id's bestimmen
            #namen sind bezogen auf das aktuelle Viereck
            bottomleft=(k-1)*nx+(k-1)+l;
            bottomright=bottomleft+1;
            topleft=k*(nx+1)+l;
            topright=topleft+1;
            middle=(nx+1)*(ny+1)+(k-1)*nx+l

            #append!(incf,[bottomleft,  bottomright, middle])
            #append!(incf,[topleft,     topright,    middle])
            #append!(incf,[bottomright, topright,    middle])
            #append!(incf,[bottomleft,  topleft,     middle])

            if mod(k,2)==1
                append!(incf,[bottomleft,  bottomright, middle])
                append!(incf,[middle,      topright,    topleft])
                append!(incf,[bottomright, topright,    middle])
                append!(incf,[middle,      topleft,     bottomleft])
            else
                append!(incf,[middle,     bottomright, bottomleft])
                append!(incf,[topleft,    topright,    middle])
                append!(incf,[middle,     topright,    bottomright])
                append!(incf,[bottomleft, topleft,     middle])
            end
        end
    end

    #Berechnen der Inzidenz 1->0
    ince=Int[];
    z=1;

    for k in 1:(nx+1)*ny
        i=[z, z+nx+1];
        append!(ince,i);
        z+=1;
    end
    zh=1;
    for h in 1:nx
        z=zh;
        for k in 1:(ny+1)
            i=[z, z+1]
            append!(ince,i);
            z+=nx+1
        end
        zh+=1;
    end

    #schräge Kanten
    #reihenfolge: (immer von ecke nach mitte) unten links, unten rechts, oben rechts, oben links
    z=1;
    for k in 1:ny
        for l in 1:nx
            bottomleft=(k-1)*nx+(k-1)+l;
            bottomright=bottomleft+1;
            topleft=k*(nx+1)+l;
            topright=topleft+1;
            middle=(nx+1)*(ny+1)+(k-1)*nx+l

            append!(ince,[bottomleft, middle])
            append!(ince,[bottomright, middle])
            append!(ince,[topleft, middle])
            append!(ince,[topright, middle])
        end
    end

    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>incf,"10"=>ince);

    #Initialisieren des Boundaryvektoren mit Einträgen
    #   0 für innere Kante/Knoten,
    #   1 für RandKante/Knoten mit freeslip,
    #   -x für periodische RandKante/Knoten mit x gegenüberliegende RandKante/Knoten
    bE=spzeros(Int, size[2]);
    bV=spzeros(Int, size[1]);
    ordV=Int[];
    corners=Set{Int}();
    #Berechnen des Boundaryvektors für die Randknoten
    for i in 1:size[1]
        if (coord[1,i]==xl || coord[1,i]==xr) && (coord[2,i]==yl || coord[2,i]==yr)
            push!(corners,i);
            if condEW==:constant && condTB==:constant
                bV[i]=1;
            elseif condEW==:periodic && condTB==:periodic && !(coord[1,i]==xl && coord[2,i]==yl)
                bV[i]=-1;
            elseif condEW==:periodic && condTB==:constant && coord[1,i]==xl
                bV[i]=-i-nx;
            elseif condTB==:periodic && condEW==:constant && coord[2,i]==yl
                bV[i]=-i-(nx+1)*ny;
            end
        elseif condEW==:constant && (coord[1,i]==xl || coord[1,i]==xr)
            bV[i]=1;
        elseif coord[1,i]==xl
            #condEW==:periodic
            bV[i]=-i-nx;
        elseif condTB==:constant && (coord[2,i]==yl || coord[2,i]==yr)
            bV[i]=1;
        elseif coord[2,i]==yl
            #condTB==:periodic
            bV[i]=-i-(nx+1)*ny;
        end
    end

    #Berechnen des Boundaryvektors für die Randkanten
    for i in 1:size[2]
        v1=ince[offe[i]];
        v2=ince[offe[i]+1]
        b1=bV[v1];
        b2=bV[v2];

        #bei schrägen Kanten continue
        if !(isequal(coord[1,v1],coord[1,v2]) || isequal(coord[2,v1],coord[2,v2]))
            continue;
        end

        if b1==0 && b2==0
            continue;
        elseif v1==1 || v2==1 #Kanten unten links gesondert behandelt wegen periodic/periodic, 0 in bV[1] macht Probleme
            if b2==b1==1
                bE[i]=1;
            elseif (b1==1 && b2<=0) || (b1<=0 && b2==1)
                bE[i]=1;
            elseif (b1==0 && b2<0) || (b1<0 && b2==0) || (b1<0 && b2<0)
                if coord[1,v1]==coord[1,v2]
                    bE[i]=-i-nx;
                elseif coord[2,v1]==coord[2,v2]
                    bE[i]=-i-ny;
                end
            end
        elseif (b1!=0 && b2!=0) || in(v1,corners) || in(v2,corners)
            if b2==b1==1
                bE[i]=1;
            elseif (b1==1 && b2<=0) || (b1<=0 && b2==1)
                bE[i]=1;
            elseif b1<0 && b2<0
                if coord[1,v1]==coord[1,v2]
                    bE[i]=-i-nx;
                elseif coord[2,v1]==coord[2,v2]
                    bE[i]=-i-ny;
                end
            end
        end
    end


    #Initialisieren der Topologie, Geometrie und damit des Meshes
    n=Int[nx,ny];
    l=Float64[xl,yl];
    r=Float64[xr,yr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV);

    return m
end



#Funktion zum Generieren von zweidimensionalen, gleichschenkligen Dreieck-Gitter
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       xl bzw. yl ist die Länge des Meshes in x- bzw. y-Richtung, als Default ist die Länge nx bzw. ny
function generateTriMesh3(nx::Int, ny::Int, xl::Float64=0.0, yl::Float64=0.0, xr::Float64=Float64(nx), yr::Float64=Float64(ny))
    if iseven(ny)
        #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
        nf=ny*(2*nx-1);
        ne=(ny/2+1)*nx+ny/2*(nx-1)+ny*2*nx
        nv=ny/2*(nx+1)+ny/2*nx+nx+1;
    else
        nf=ny*(2*nx-1);
        ne=(ny+1)/2*nx+(ny+1)/2*(nx-1)+ny*2*nx;
        nv=(ny+1)/2*(nx+1)+(ny+1)/2*nx;
    end

    size=[Int(nv), Int(ne), Int(nf)];
    nk=3;

    #Initialisieren des Offsets mit den Einträgen "20" und "10"
    off=Dict("20"=>collect(1:nk:(nk*size[3]+1)),"10"=>collect(1:2:(2*size[2]+1)));

    #Berechnen der Inzidenz 2->0
    incf=Int[];
    z=1;
    for k in 1:ny
        iseven(k) ? e=nx-1 : e=nx
        if isodd(k)
            for h in 1:(nx-1)
                i=[z, z+1, z+1+nx];
                append!(incf,i);
                i=[z+1, z+2+nx, z+1+nx];
                append!(incf,i);
                z+=1;
            end
            i=[z, z+1, z+1+nx];
            append!(incf,i);
            z+=2;
        else
            for h in 1:(nx-1)
                i=[z, z+1+nx, z+nx];
                append!(incf,i);
                i=[z, z+1, z+1+nx];
                append!(incf,i);
                z+=1;
            end
            i=[z, z+1+nx, z+nx];
            append!(incf,i);
            z+=1;
        end

    end


    #Berechnen der Inzidenz 1->0
    ince=Int[];
    z=1;
    for k in 1:(ny+1)
        isodd(k) ? e=nx : e=nx-1;
        for h in 1:e
            i=[z, z+1];
            append!(ince,i);
            z+=1;
        end
        z+=1;
    end
    z=1;

    for k in 1:ny
        if isodd(k)
            for h in 1:nx
                i=[z, z+nx+1];
                append!(ince,i);
                i=[z+1, z+nx+1];
                append!(ince,i);
                z+=1;
            end
            z+=1
        else
            for h in 1:nx
                i=[z, z+nx];
                append!(ince,i);
                i=[z, z+nx+1];
                append!(ince,i);
                z+=1;
            end
        end
    end

    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>incf,"10"=>ince);

    #Berechnen der Koordinatenmatrix, basierend auf einem äquidistanten Gitter
    coord=Array{Float64}(undef,2,size[1]);
    ax=(xr-xl)/nx;
    ay=(yr-yl)/ny;
    zx=xl;
    zy=yl;
    z=1;
    for k in 1:(ny+1)
        if isodd(k)
            for l in 1:(nx+1)
                coord[:,z]=[zx, zy];
                zx+=ax;
                z+=1;
            end
            zx=xl+ax/2;
        else
            for l in 1:nx
                coord[:,z]=[zx, zy];
                zx+=ax;
                z+=1;
            end
            zx=xl;
        end
        zy+=ay;
    end

    #Initialisieren der Topologie, Geometrie und damit des Meshes
    n=Int[nx,ny];
    l=Float64[xl,yl];
    r=Float64[xr,yr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG);

    return m
end


#Funktion zum Generieren von zweidimensionalen Sechseck-Gitter
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       Die Länge jeder Seite jedes Sechseck ist hier konstant 1.
function generateHexMesh(nx::Int, ny::Int)
    if iseven(ny)
        #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
        nf=0.5*ny*(2*nx-1)
        ne=3*nx*ny+0.5*ny+2*(nx-1) # =2*nx*ny+0.5*ny*(nx+1)+0.5*ny*nx+nx
        nv=ny*(2*nx+1)+2*nx-1

        size=[Int(nv), Int(ne), Int(nf)];
        nk=6;

        #Initialisieren des Offsets mit den Einträgen "20" und "10"
        off=Dict("20"=>Array(1:nk:(nk*size[3]+1)),"10"=>Array(1:2:(2*size[2]+1)));

        #Berechnen der Inzidenz 2->0
        incf=Int[];
        z=1;
        for k in 1:(ny-1)
            iseven(k) ? e=nx-1 : e=nx
            for h in 1:e
                i=[z, z+1, z+2, z+2*nx+3, z+2*nx+2, z+2*nx+1];
                append!(incf,i);
                z+=2;
            end
            z+=2;
        end

        for h in 1:(nx-1)
            i=[z, z+1, z+2, z+2*nx+2, z+2*nx+1, z+2*nx];
            append!(incf,i);
            z+=2;
        end

        #Berechnen der Inzidenz 1->0
        ince=Int[];
        z=1;
        for k in 1:ny
            for h in 1:2*nx
                i=[z, z+1];
                append!(ince,i);
                z+=1;
            end
            z+=1;
        end
        for h in 1:2*(nx-1)
            i=[z, z+1];
            append!(ince,i);
            z+=1;
        end

        z=1;
        for k in 1:(nx+1)*ny*0.5+0.5*ny*nx-nx
            i=[z, z+2*nx+1];
            append!(ince,i);
            z+=2;
        end
        for k in 1:nx
            i=[z, z+2*nx];
            append!(ince,i);
            z+=2;
        end

        #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
        inc=Dict("20"=>incf,"10"=>ince);

        #Berechnen der Koordinatenmatrix, basierend auf einem äquidistanten Gitter
        coord=Array{Float64}(undef,2,size[1]);
        ax=1;
        ay=1;
        zx=0.0;
        zy=ay;

        z=1;
        for k in 1:ny
            for l in 1:(2*nx+1)
                coord[:,z]=[zx, zy];
                zx+=ax;
                if isodd(k)
                    iseven(l) ? zy+=ay : zy-=ay;
                else
                    iseven(l) ? zy-=ay : zy+=ay;
                end
                z+=1;
            end
            zx=0.0;
            zy+=2*ay;
        end

        zx=ax;
        zy-=ay;
        for k in 1:(2*nx-1)
            coord[:,z]=[zx, zy];
            zx+=ax;
            iseven(k) ? zy-=ay : zy+=ay;
            z+=1;
        end

        #Initialisieren der Topologie, Geometrie und damit des Meshes
        n=Int[nx,ny];
        l=Float64[0.0,0.0];
        r=Float64[Float64(nx),Float64(ny)];
        mT=meshTopology(inc,off,n);
        mG=meshGeometry(coord,l,r);
        m=mesh(mT,mG);

        return m
    else
        #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
        nf=0.5*(ny+1)*nx+0.5*(ny-1)*(nx-1)
        ne=2*nx*(ny+1)+(ny-1)*0.5*nx+(ny+1)*0.5*(nx+1)
        nv=(ny+1)*(2*nx+1)

        size=[Int(nv), Int(ne), Int(nf)];
        nk=6;

        #Initialisieren des Offsets mit den Einträgen "20" und "10"
        off=Dict("20"=>Array(1:nk:(nk*size[3]+1)),"10"=>Array(1:2:(2*size[2]+1)));

        #Berechnen der Inzidenz 2->0
        incf=Int[];
        z=1;
        for k in 1:ny
            iseven(k) ? e=nx-1 : e=nx
            for h in 1:e
                i=[z, z+1, z+2, z+2*nx+3, z+2*nx+2, z+2*nx+1];
                append!(incf,i);
                z+=2;
            end
            z+=2;
        end

        #Berechnen der Inzidenz 1->0
        ince=Int[];
        z=1;
        for k in 1:(ny+1)
            for h in 1:2*nx
                i=[z, z+1];
                append!(ince,i);
                z+=1;
            end
            z+=1;
        end
        z=1;
        for k in 1:(nx+1)*(ny+1)*0.5+0.5*(ny-1)*nx
            i=[z, z+2*nx+1];
            append!(ince,i);
            z+=2;
        end

        #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
        inc=Dict("20"=>incf,"10"=>ince);

        #Berechnen der Koordinatenmatrix, basierend auf einem äquidistanten Gitter
        coord=Array{Float64}(undef,2,size[1]);
        ax=1;
        ay=1;
        zx=0.0;
        zy=ay;
        z=1;
        for k in 1:(ny+1)
            for l in 1:(2*nx+1)
                coord[:,z]=[zx, zy];
                zx+=ax;
                if isodd(k)
                    iseven(l) ? zy+=ay : zy-=ay;
                else
                    iseven(l) ? zy-=ay : zy+=ay;
                end
                z+=1;
            end
            zx=0.0;
            zy+=2*ay;
        end

        #Initialisieren der Topologie, Geometrie und damit des Meshes
        n=Int[nx,ny];
        l=Float64[0.0,0.0];
        r=Float64[Float64(nx),Float64(ny)];
        mT=meshTopology(inc,off,n);
        mG=meshGeometry(coord,l,r);
        m=mesh(mT,mG);

        return m
    end
end

function generateCubMesh(nx::Int, ny::Int, nz::Int, condEW::Symbol, condTB::Symbol, xl::Float64=0.0, yl::Float64=0.0, zl::Float64=0.0, xr::Float64=Float64(nx), yr::Float64=Float64(ny), zr::Float64=Float64(nz))

    #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
    size=[(nz+1)*(ny+1)*(nx+1), (ny*(nx+1)+nx*(ny+1))*(nz+1)+nz*(nx+1)*(ny+1), (nz+1)*nx*ny+nz*(ny*(nx+1)+nx*(ny+1)), nz*nx*ny];
    nk=4;

    #Initialisieren des Offsets mit den Einträgen "20", "10" und "30"
    off=Dict("30"=>collect(1:8:(8*size[4]+1)),"20"=>collect(1:4:(4*size[3]+1)),"10"=>collect(1:2:(2*size[2]+1)));

    #Berechnen der Koordinatenmatrix, basierend auf einem äquidistanten Gitter
    coord=Array{Float64,2}(undef,3,size[1]);
    ax=(xr-xl)/nx;
    ay=(yr-yl)/ny;
    az=(zr-zl)/nz;
    zx=xl;
    zy=yl;
    zz=zl;
    z=1;
    for j in 1:(nz+1)
        for k in 1:(ny+1)
            for l in 1:(nx+1)
                coord[:,z]=[zx, zy, zz];
                zx+=ax;
                z+=1;
            end
            zx=xl;
            zy+=ay;
        end
        zy=yl;
        zz+=az;
    end

    #Berechnen der Inzidenz 3->0
    incc=Int[];
    hz=(nx+1)*(ny+1);
    z=1;
    for j in 1:nz
        for k in 1:ny
            for h in 1:nx
                i=[z, z+1, z+1+nx+1, z+nx+1, z+hz, z+hz+1, z+hz+nx+2, z+hz+nx+1];
                append!(incc,i);
                z+=1;
            end
            z+=1;
        end
        z+=nx+1;
    end

    #Berechnen der Inzidenz 2->0
    incf=Int[];
    z=1;
    for j in 1:(nz+1)
        for k in 1:ny
            for h in 1:nx
                i=[z, z+1, z+1+nx+1, z+nx+1];
                append!(incf,i);
                z+=1;
            end
            z+=1;
        end
        z+=nx+1;
    end
    z=1;
    for j in 1:nz
        wz=copy(z);
        for k in 1:(ny+1)
            for h in 1:nx
                i=[z, z+1, z+hz+1, z+hz];
                append!(incf,i);
                z+=1;
            end
            z+=1;
        end
        z=wz;
        for k in 1:ny
            for h in 1:(nx+1)
                i=[z, z+nx+1, z+hz+nx+1, z+hz];
                append!(incf,i);
                z+=1;
            end
            #z+=1;
        end
        z+=nx+1;
    end

    #Berechnen der Inzidenz 1->0
    ince=Int[];
    z=1;
    for j in 1:(nz+1)
        wz=copy(z);
        for k in 1:(ny+1)
            for h in 1:nx
                i=[z, z+1];
                append!(ince,i);
                z+=1;
            end
            z+=1;
        end
        z=wz;
        for k in 1:(nx+1)*ny
            i=[z, z+nx+1];
            append!(ince,i);
            z+=1;
        end
        z+=nx+1;
    end
    z=1;
    for j in 1:nz
        for k in 1:(ny+1)
            for h in 1:(nx+1)
                i=[z, z+hz];
                append!(ince,i);
                z+=1;
            end
        end
    end


    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("30"=>incc,"20"=>incf,"10"=>ince);

    #Initialisieren des Boundaryvektoren mit Einträgen
    #   0 für innere Kante/Knoten,
    #   1 für RandKante/Knoten mit freeslip,
    #   -x für periodische RandKante/Knoten mit x gegenüberliegende RandKante/Knoten
    bE=spzeros(Int, size[2]);
    bV=spzeros(Int, size[1]);
    #Randbetrachtung fehlt

    #Initialisieren der Topologie, Geometrie und damit des Meshes
    n=Int[nx,ny,nz];
    l=Float64[xl,yl,zl];
    r=Float64[xr,yr,zr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV);

    return m
end
