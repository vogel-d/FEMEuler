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
  D::Int
  n::Array{Int,1}
end

function meshTopology(inc::Dict{String, Array{Int,1}},off::Dict{String, Array{Int,1}}, n::Array{Int,1})
  size=Int[];
  push!(size,maximum(inc["10"]));
  a=0;
  for k in 2:(length(inc)+1)
    push!(size,length(off["$(k-1)$a"])-1);
  end
  D=length(size)-1;
  meshTopology(inc,off,size,D,n)
end


#Struct-Objekt, dass die Geometrie eines Meshes speichert durch das
#Speichern der Koordinaten der vertices in einer Matrix, in der in jeder Spalte i
#die Koordinaten für den vertex i gespeichert sind
#l und r speichern jeweils die äußersten Koordinaten für jede Dimension
struct meshGeometry
  coordinates::Array{AbstractFloat,2}
  l::Array{AbstractFloat,1}
  r::Array{AbstractFloat,1}
end

#Struct-Objekt zum Speichern eines Meshes, durch Speichern der
#Topologie als meshTopology-Objekt sowie der Geometrie als meshGeometry-Objekt
struct mesh
  topology::meshTopology
  geometry::meshGeometry
  meshType::Int
  edgeLength::Array{AbstractFloat,1} #Kantenlängen
  normals::Array{AbstractFloat,2} #Normalen des Referenzelementes
  boundaryEdges::SparseVector{Int,Int};
  boundaryVertices::SparseVector{Int,Int};
end

function mesh(topology::meshTopology, geometry::meshGeometry, bE::SparseVector{Int,Int}, bV::SparseVector{Int,Int})
  inc=topology.incidence["10"];
  coord=geometry.coordinates;
  ne=topology.size[2];
  l=Array{AbstractFloat,1}(undef,ne);
  z=1;
  for i in 1:ne
      c=coord[:,inc[z:(z+1)]];
      l[i]=sqrt((c[1,1]-c[1,2])^2+(c[2,1]-c[2,2])^2);
      z+=2;
  end
  mt=topology.offset["20"][2]-topology.offset["20"][1];
  mt==4 ? n=[0.0 -1.0 0.0 1.0;-1.0 0.0 1.0 0.0] : n=[0.0 -1.0 0.7071067811865475244;-1.0 0.0 0.7071067811865475244];
  mesh(topology, geometry, mt, l, n, bE, bV)
end
