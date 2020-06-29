#Funktion zum Generieren von zweidimensionalen Rechteck-Gittern
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       xl bzw. yl ist die Länge des Meshes in x- bzw. y-Richtung, als Default ist die Länge nx bzw. ny
function generateRectMesh(nx::Int, ny::Int, condEW::Symbol, condTB::Symbol, xl::Float64=0.0, xr::Float64=Float64(nx), yl::Float64=0.0, yr::Float64=Float64(ny); meshType::Int64=4)

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
    m=mesh(mT,mG, bE, bV,condEW,condTB,meshType);

    return m
end

#Funktion zum Generieren von zweidimensionalen Dreieck-Gitter im Rahmen eines Rechtecks
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       xl bzw. yl ist die Länge des Meshes in x- bzw. y-Richtung, als Default ist die Länge nx bzw. ny
function generateTriMeshHalved(nx::Int, ny::Int, condEW::Symbol, condTB::Symbol, xl::Float64=0.0, xr::Float64=Float64(nx), yl::Float64=0.0, yr::Float64=Float64(ny))

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
    m=mesh(mT,mG, bE, bV,condEW,condTB);

    return m
end

function generateTriMeshQuartered(nx::Int, ny::Int, condEW::Symbol, condTB::Symbol, xl::Float64=0.0, xr::Float64=Float64(nx), yl::Float64=0.0, yr::Float64=Float64(ny))

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
            #variablennamen sind bezogen auf das aktuelle Viereck
            bottomleft=(k-1)*nx+(k-1)+l;
            bottomright=bottomleft+1;
            topleft=k*(nx+1)+l;
            topright=topleft+1;
            middle=(nx+1)*(ny+1)+(k-1)*nx+l;

            #nummeriert nach aufsteigend globaler Nummer
            append!(incf,[bottomleft,  bottomright, middle])
            append!(incf,[bottomright, topright,    middle])
            append!(incf,[topleft,     topright,    middle])
            append!(incf,[bottomleft,  topleft,     middle])
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

            append!(ince,[bottomleft,  middle])
            append!(ince,[bottomright, middle])
            append!(ince,[topright,    middle])
            append!(ince,[topleft,     middle])
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
    m=mesh(mT,mG, bE, bV,condEW,condTB);

    return m
end

#Funktion zum Generieren von zweidimensionalen, gleichschenkligen Dreieck-Gitter
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       xl bzw. yl ist die Länge des Meshes in x- bzw. y-Richtung, als Default ist die Länge nx bzw. ny
function generateTriMeshEquilateral(xl::Float64, xr::Float64, yl::Float64, yr::Float64, nrows::Int64, condEW::Symbol, condTB::Symbol)
    l = (2*(yr-yl))/(nrows*sqrt(3));
    height=(yr-yl)/nrows;
    nx=Int64(div(xr-xl-0.5*l,l)+1); #smallest nx so that nx*l+0.5*l>yr
    xR=xl+nx*l+0.5*l;
    @info "mesh details\n ny=$nrows (=nrows)\n nx=$nx \n edge length=$l\n gridsize=[$xl,$xR]x[$yl,$yr]\n xmax deviation: $(xR-xr)"

    #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
    size=[(nx+1)*(nrows+1), nrows*(3*nx+1)+nx, nrows*2*nx];
    nk=3;

    #Initialisieren des Offsets mit den Einträgen "20" und "10"
    offe=collect(1:2:(2*size[2]+1));
    off=Dict("20"=>collect(1:nk:(nk*size[3]+1)),"10"=>offe);

    #Erstellen der Koordinatenmatrix
    coord=Array{Float64}(undef,2,size[1]);
    for line in 1:(nrows+1)
        #coordinates for one horizontal line of points are assembled at once
        coord[1,((line-1)*(nx+1)+1):((line-1)*(nx+1)+1+nx)] = collect(range(mod(line-1,2)*(0.5*l)+xl,step=l,length=nx+1))
        coord[2,((line-1)*(nx+1)+1):((line-1)*(nx+1)+1+nx)] = ones(nx+1)*(yl+(line-1)*height)
    end

    #note to mod(,2) if-clauses: the two types of rows differ in calculating
    #                            the vertex numbers of vertices above

    #Erstellen der Inzidenz 2->0
    inc20=Int[];
    for row in 1:nrows
        if mod(row,2)==1
            for vertex in range((nx+1)*(row-1)+1, step=1, length=nx)
                append!(inc20,[vertex, vertex+1, vertex+nx+1])
                append!(inc20,[vertex+1, vertex+1+nx, vertex+1+nx+1])
            end
        else
            for vertex in range((nx+1)*(row-1)+1, step=1, length=nx)
                append!(inc20,[vertex, vertex+nx+1, vertex+nx+2])
                append!(inc20,[vertex, vertex+1, vertex+1+nx+1])
            end
        end
    end

    #Erstellen der Inzidenz 1->0
    inc10=Int[];
    for vertex in 1:((nx+1)*(nrows+1))
        mod(vertex,nx+1)==0 && continue; #skip the last vertex of each horizontal line
        append!(inc10,[vertex,vertex+1]);
    end
    for row in 1:nrows
        if mod(row,2)==1
            bottomleft=(nx+1)*(row-1)+1;
            bottomright=bottomleft+1;
            top=bottomright+nx;
            for hat_number in 0:(nx-1)
                append!(inc10,[bottomleft,top]);
                append!(inc10,[bottomright,top]);
                bottomleft+=1;
                bottomright+=1;
                top+=1;
            end
            #add last edge of row of odd number
            append!(inc10,[bottomleft,top]);
        else
            bottom=(nx+1)*(row-1)+1;
            topleft=bottom+nx+1;
            topright=topleft+1;
            for hat_number in 0:(nx-1)
                append!(inc10,[bottom,topleft]);
                append!(inc10,[bottom,topright]);
                bottom+=1;
                topleft+=1;
                topright+=1;
            end
            #add last edge of row of even number
            append!(inc10,[bottom,topleft]);
        end
    end

    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>inc20,"10"=>inc10);


    bE=spzeros(Int, size[2]);
    bV=spzeros(Int, size[1]);

    #boundary-Vektor Randknoten
    if condTB==:constant
        bV[1:nx+1]=ones(nx+1);
        bV[(size[1]-nx):size[1]]=ones(nx+1);
    elseif condTB==:periodic
        bV[1:nx+1]=collect(-(size[1]-nx):-1:-size[1]);
    end

    if condEW==:constant
        bV[1:(nx+1):(size[1]-nx)]=ones(nrows+1);
        bV[(nx+1):(nx+1):size[1]]=ones(nrows+1);
    elseif condEW==:periodic
        bV[1:(nx+1):(size[1]-nx)]=-(nx+1):-(nx+1):-size[1];
    end


    #boundary-Vektor Randkanten
    if condTB==:constant
        bE[1:nx]=ones(nx);
        bE[(nrows*nx+1):((nrows+1)*nx)]=ones(nx);
    elseif condTB==:periodic
        bE[1:nx]=-(nrows*nx+1):-1:-((nrows+1)*nx);
    end

    if condEW==:constant
        skipHorizontal=(nrows+1)*nx;
        bE[range(skipHorizontal+1,step=2*nx+1,length=nrows)]=ones(nrows);
        bE[range(skipHorizontal+1+2*nx,step=2*nx+1,length=nrows)]=ones(nrows);
    elseif condEW==:periodic
        skipHorizontal=(nrows+1)*nx;
        bE[range(skipHorizontal+1,step=2*nx+1,length=nrows)]=range(-(skipHorizontal+1+2*nx),step=-(2*nx+1),length=nrows);
    end

    n=Int[nx,nrows];
    l=Float64[xl,yl];
    r=Float64[xR,yr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV,condEW,condTB);

    return m
end


#Funktion zum Generieren von zweidimensionalen, gleichschenkligen Dreieck-Gitter
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       xl bzw. yl ist die Länge des Meshes in x- bzw. y-Richtung, als Default ist die Länge nx bzw. ny
function generateTriMeshIsosceles(nx::Int, ny::Int, xl::Float64=0.0, yl::Float64=0.0, xr::Float64=Float64(nx), yr::Float64=Float64(ny))
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

function generateHexMesh(xl::Float64, xr::Float64, yl::Float64, yr::Float64, nrows::Int64, condEW::Symbol, condTB::Symbol; meshType::Int64=4)
    (isodd(nrows) && condTB==:periodic) && error("vertical periodic boundary only possible for even nrows.")

    l = ((yr-yl)/nrows) * (2/3);

    yR=yl+(3/2)*nrows*l+0.5*l;

    hx=sqrt(3)*l/2
    nx=Int(ceil((xr-xl-hx)/(2*hx)));
    xR=xl+(2*nx+1)*hx;

    #Berechnen der Anzahl der Entitäten für die verschiedenen Dimensionen
    size=Int64[(nrows+1)*(2*(nx+1))-2, nrows*(nx+1)+(nrows+1)*(2*nx+1)-2, nx*nrows];
    nk=6;

    @info "mesh details\n ny=$nrows (=nrows)\n nx=$nx \n edge length=$l\n gridsize=[$xl,$xR]x[$yl,$yR]\n xmax deviation: $(xR-xr)\n ymax deviation: $(yR-yr)\n mesh type=$meshType for use of compound elements"


    #Initialisieren des Offsets mit den Einträgen "20" und "10"
    offe=collect(1:2:(2*size[2]+1));
    off=Dict("20"=>collect(1:nk:(nk*size[3]+1)),"10"=>offe);

    #Erstellen der Koordinatenmatrix
    coord=Array{Float64}(undef,2,size[1]);
    coordx=Float64[];
    coordy=Float64[];
    line_ypos=-l;
    sign=1.0;

    for line in 1:(nrows+1)
        #get line parameters
        line_ypos+=l+iseven(line)*l;
            #(hats pointing up or down)
        isodd(line) ? sign=1.0 : sign=-1.0

        #set first vertex of line (except first line and last line if nrows is odd)
        if line!=1 && ((line!=nrows+1) || iseven(nrows))
            append!(coordx,xl)
            append!(coordy,line_ypos)
        end

        #set second vertex of line (beginning of first hat)
        append!(coordx,xl+hx)
        append!(coordy,line_ypos+sign*0.5*l)

        #set nx hats (upside down depending on even or odd row)
        for hat in 1:nx
            append!(coordx,xl+2*hat*hx)
            append!(coordy,line_ypos)

            append!(coordx,xl+2*hat*hx+hx)
            append!(coordy,line_ypos+sign*0.5*l)
        end
    end
    #delete the very last vertex if nrows is even,
    #did so for hat-constructing-continuity
    if iseven(nrows)
        pop!(coordx)
        pop!(coordy)
    end
    coord[1,:]=coordx;
    coord[2,:]=coordy;

    #Erstellen der Inzidenz 2->0
    inc20=Int[];
    bottomleft=-2
    for row in 1:nrows
        bottomleft+=1+2*isodd(row)
        for cell in 1:nx
            #topleft is usually bottomleft+2nx+2, only the last row will miss one vertex if odd
            topleft=bottomleft+2*nx+1+1*(row!=nrows || iseven(nrows))
            append!(inc20,[bottomleft,bottomleft+1,bottomleft+2,topleft+2,topleft+1,topleft])
            bottomleft+=2
        end
    end

    #Erstellen der Inzidenz 1->0
    inc10=Int[]
    #vertical edges
    bottom=-2;
    top=0;
    for row in 1:nrows
        bottom+=1+2*isodd(row)
        for edge in 1:nx
            top=bottom+2*nx+1+1*(row!=nrows || iseven(nrows))
            append!(inc10,[bottom,top])
            bottom+=2
        end
        top=bottom+2*nx+1+1*(row!=nrows || iseven(nrows))
        append!(inc10,[bottom,top])
    end
    #rest of edges (every vertex has to be connected to its number neighbors except last in lines)
    for vertex in 1:size[1]-1
        #each line has 2*nx vertices for the hats plus 2 outer vertices
        #the very first left outer vertex is missing, so its mod(vertex PLUS ONE,2nx+2) (missing vertex results in right outer vertices having 2*nx+2-1)
        if mod(vertex+1,2*nx+2)!=0
            append!(inc10,[vertex,vertex+1]);
        end
    end

    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>inc20,"10"=>inc10);


    bE=spzeros(Int, size[2]);
    bV=spzeros(Int, size[1]);

    #boundary Vektor Randknoten
    if condTB==:constant
        bV[2:(2*nx)].=1.0;
        bV[(size[1]-2*nx+1):(size[1]-1)].=1.0;
    elseif condTB==:periodic
        bV[2:(2*nx)]=-(size[1]-2*nx+2):-1:-(size[1]);
    end

    if condEW==:constant
        #walk up east/west boundary vertex-wise
        vertex_west=0;
        vertex_east=2*nx;
        for row in 1:floor(nrows/2)
            vertex_east+=1;
            vertex_west+=1;
            bV[vertex_west]=1.0;
            bV[vertex_east]=1.0;
            vertex_east+=2*nx+2;
            vertex_west+=2*nx+2;
            bV[vertex_west]=1.0;
            bV[vertex_east]=1.0;
            vertex_east+=-1;
            vertex_west+=-1;
            bV[vertex_west]=1.0;
            bV[vertex_east]=1.0;
            vertex_east+=2*nx+2;
            vertex_west+=2*nx+2;
            bV[vertex_west]=1.0;
            bV[vertex_east]=1.0;
        end
        if isodd(nrows)
            vertex_east+=1;
            vertex_west+=1;
            bV[vertex_west]=1.0;
            bV[vertex_east]=1.0;
            vertex_east+=2*nx+1;
            vertex_west+=2*nx+1;
            bV[vertex_west]=1.0;
            bV[vertex_east]=1.0;
        end
    elseif condEW==:periodic
        #walk up east/west boundary vertex-wise
        vertex_west=0;
        vertex_east=2*nx;
        for row in 1:floor(nrows/2)
            vertex_east+=1;
            vertex_west+=1;
            bV[vertex_west]=-vertex_east;
            vertex_east+=2*nx+2;
            vertex_west+=2*nx+2;
            bV[vertex_west]=-vertex_east;
            vertex_east+=-1;
            vertex_west+=-1;
            bV[vertex_west]=-vertex_east;
            vertex_east+=2*nx+2;
            vertex_west+=2*nx+2;
            bV[vertex_west]=-vertex_east;
        end
        if isodd(nrows)
            vertex_east+=1;
            vertex_west+=1;
            bV[vertex_west]=-vertex_east;
            vertex_east+=2*nx+1;
            vertex_west+=2*nx+1;
            bV[vertex_west]=-vertex_east;
        end
    end

    if condTB==:periodic && condEW==:periodic
        #glue corners together
        bV[1]=0;
        bV[2*nx+1]=-1;
        bV[size[1]-2*nx+1]=-1;
        bV[2*nx]=-size[1];
        bV[size[1]-2*nx]=-size[1];
    end

    #boundary Vektor für die Randkanten
    nVertical=nrows*(nx+1);
    if condTB==:constant
        bE[(nVertical+1):(nVertical+2*nx)].=1.0;
        bE[(size[2]-2*nx+1):size[2]].=1.0;
    elseif condTB==:periodic
        bE[nVertical+2*nx]=1.0;
        bE[size[2]-2*nx+1]=1.0;
        bE[(nVertical+1):(nVertical+2*nx-1)]=(-(size[2]-2*nx+2)):(-1):(-size[2]);
    end

    if condEW==:constant
        bE[1:(nx+1):(nVertical-nx)].=1.0
        bE[(1+nx):(nx+1):nVertical].=1.0
        bE[(nVertical+(2*nx+1)):(2*nx+1):(size[2]-4*nx)].=1.0
        bE[(nVertical+(2*nx+1)+2*nx):(2*nx+1):(size[2]-4*nx+2*nx)].=1.0
    elseif condEW==:periodic
        bE[1:(nx+1):(nVertical-nx)]=-(1+nx):-(nx+1):-nVertical
        bE[(nVertical+(2*nx+1)):(2*nx+1):(size[2]-4*nx)]=-(nVertical+(2*nx+1)+2*nx):-(2*nx+1):-(size[2]-4*nx+2*nx)
    end

    if condTB==:periodic && condEW==:periodic
        #glue corners together
        bE[size[2]-2*nx+1]=-(nVertical+2*nx)
    end

    n=Int[nx,nrows];
    l=Float64[xl,yl];
    r=Float64[xR,yR];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV,condEW,condTB, meshType);
    return m
end


#Funktion zum Generieren von zweidimensionalen Sechseck-Gitter
#Input: nx bzw. ny ist die Anzahl der Gitterelemente in x- bzw. y-Richtung, also die Feinheit des Meshes
#       Die Länge jeder Seite jedes Sechseck ist hier konstant 1.
function generateHexMesh2(nx::Int, ny::Int)
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
        bE=spzeros(Int, size[2]);
        bV=spzeros(Int, size[1]);
        m=mesh(mT,mG,bE,bV,condEW,condTB);

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
    m=mesh(mT,mG, bE, bV,condEW,condTB);

    return m
end
