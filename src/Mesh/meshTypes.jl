#Struct-Objekt, dass die Topologie eines Meshs speichert,dabei gilt:

#incidence speichert unter dem Key "ds" die Inzidenzen d->s und
#das zugehörige Offset mit dem gleichen Key speichert das Offset von incidence,
#wobei immer der Anfangsindex für jede Entität der Dimension d gespeichert wird
#und als letztes Element noch die Länge von incidence+1

#n ist die Anzahl der Zellen in jede Richtung (jede Dimension)

#size speichert für jede Dimension die Anzahl der Entitäten und D die größte Dimension
#Beide werden allerdings nicht im Konstruktor angegeben, da sie aus incidence und
#offset berechnet werden
struct meshTopology
  incidence::Dict{String, Array{Int,1}}
  offset::Dict{String, Array{Int,1}}
  size::Array{Int,1}
  dim::Int
  n::Array{Int,1}
end

function meshTopology(inc::Dict{String, Array{Int,1}},off::Dict{String, Array{Int,1}}, n::Array{Int,1})
  size=Int[];
  push!(size,maximum(inc["10"]));
  a=0;
  for k in 2:(length(inc)+1)
    push!(size,length(off["$(k-1)$a"])-1);
  end
  dim=length(size)-1;
  meshTopology(inc,off,size,dim,n)
end


#Struct-Objekt, dass die Geometrie eines Meshes speichert durch das
#Speichern der Koordinaten der vertices in einer Matrix, in der in jeder Spalte i
#die Koordinaten für den vertex i gespeichert sind
#l und r speichern jeweils die äußersten Koordinaten für jede Dimension
struct meshGeometry
  coordinates::Array{Float64,2}
  dim::Int;
  l::Array{Float64,1}
  r::Array{Float64,1}
end

function meshGeometry(coord::Array{Float64,2},l::Array{Float64,1}, r::Array{Float64,1})
  dim=size(coord,1);
  meshGeometry(coord,dim,l,r)
end

#Struct-Objekt zum Speichern eines Meshes, durch Speichern der
#Topologie als meshTopology-Objekt sowie der Geometrie als meshGeometry-Objekt
struct mesh
  topology::meshTopology
  geometry::meshGeometry
  meshType::Int
  edgeLength::Array{Float64,1} #Kantenlängen
  normals::Array{Float64,2} #Normalen des Referenzelementes
  boundaryEdges::SparseVector{Int,Int};
  boundaryVertices::SparseVector{Int,Int};
  orientation::Array{Float64,1};
  boundaryConditionEW::Symbol;
  boundaryConditionTB::Symbol;
end

function mesh(topology::meshTopology, geometry::meshGeometry, bE::SparseVector{Int,Int}, bV::SparseVector{Int,Int},condEW::Symbol,condTB::Symbol, mt::Int=topology.offset["20"][2]-topology.offset["20"][1])
  inc=topology.incidence["10"];
  coord=geometry.coordinates;
  ne=topology.size[2];
  l=Array{Float64,1}(undef,ne);
  z=1;
  for i in 1:ne
      c=coord[:,inc[z:(z+1)]];
      l[i]=sqrt(sum((c[:,1].-c[:,2]).^2));
      z+=2;
  end
  if mt==4
    if geometry.dim==3
      n=[0.0 1.0 0.0 1.0;-1.0 0.0 -1.0 0.0]
    else
      #n=[0.0 1.0 0.0 1.0;-1.0 0.0 -1.0 0.0]
      n=[0.0 1.0 0.0 1.0;1.0 0.0 1.0 0.0]
    end
  else
    n=[0.0 1.0 1.0;-1.0 1.0 0.0];
  end
  orientation=Float64[];
  mesh(topology, geometry, mt, l, n, bE, bV, orientation, condEW, condTB)
end
