#Funktion zum Berechnen der Inzidenz d->d' aus d'->d mit d'>d
function meshTranspose(m::mesh, d::Int, ds::Int)
  if d<ds
    #Bestimmen der Anzahl an Entitäten von der Dimension d bzw. d'
    nd=m.topology.size[d+1];
    nds=m.topology.size[ds+1];

    #Bestimmen der Inzidenz und des Offsets von d'->d
    ind=m.topology.incidence["$ds$d"];
    off=m.topology.offset["$ds$d"];

    z=1;
    o=Int[1];
    i=Int[];

    for k in 1:nd
      for h in 1:nds
        #Überprüfen, ob die Entität k der Dimension d inzident ist
        #zur Entität h der Dimension d', falls ja dann wird h für k zu d->d' hinzugefügt
        if any(isequal(k), ind[off[h]:off[h+1]-1])
          push!(i,h);
          z+=1;
        end
      end

      #Erweitern des offsets von d->d' mittels des Zählers z, der immer hochzählt,
      #falls ein h zu d->d' hinzugefügt wird
      push!(o,z);
    end

    return i,o
  else
    error("Zum Berechnen der Inzidenz d->d' durch meshTranspose muss gelten d<d'!");
  end
end

#Funktion zum Berchnen von d->d' aus d->d'' und d''->d' mit d>=d'
function meshIntersection(m::mesh,d::Int,ds::Int,dss::Int)
  if d>=ds
    #Berechnen der Anzahl der Entitäten der Dim. d sowie der Inzidenz und Offsets
    #von d->d'' und d''->d
    nd=m.topology.size[d+1];

    ind1=m.topology.incidence["$d$dss"];
    off1=m.topology.offset["$d$dss"];

    ind2=m.topology.incidence["$dss$ds"];
    off2=m.topology.offset["$dss$ds"];

    z=1;
    t=false;
    o=Int[1];
    i=Int[];

    #Falls die Dim d und d' gleich sind, müssen nicht die gemeinsamen vertices der
    #Entitäten überprüft werden, weshalb dies gesondert betrachtet wird
    #Insbesondere muss dies gesondert betrachtet werden, da beim Betrachten der vertices
    #für d=d'=0 auf 0->0 zugegriffen werden müsste, was nicht Teil der Eingabe ist
    if d==ds
      for k in 1:nd
        done=Int[];
        interim=Int[];
        #Bestimmen der Entitäten h der Dim d'', die inzident zu k der Dim d sind,
        #sowie der Entitäten j der Dim d', die inzident zu h sind
        for h in ind1[off1[k]:off1[k+1]-1]
          for j in ind2[off2[h]:off2[h+1]-1]
            #Da ein j inzident zu versch. h sein kann, wird durch einen Array, der die
            #bereits hinzugefügten j speichert, abgefangen dass ein j in der Inzidenz der
            #Entität k mehrmals vorkommen kann
            if in(j,done)
              continue;
            end
            #Abfangen das eine Entität inzident zu sich selbst ist
            if k!=j
              push!(interim,j);
              z+=1;
              t=true;
              push!(done,j);
            end
          end
        end
      sort!(interim);
      append!(i,interim);
      #Durch t wird verhindert, dass der Zähler,
      #der sich nicht verändert hat, zum offset hinzugefügt wird
      if t
        push!(o,z);
      end
      t=false;
      end
    #Vorgehen analog zum ersten Fall, nur wird hier zusätzlich d->0 und d'->0 bestimmt
    #und es wird überprüft, ob d'->0 eine Teilmenge von d->0 ist
    else
      a=0;
      ind3=m.topology.incidence["$d$a"]
      off3=m.topology.offset["$d$a"]
      ind4=m.topology.incidence["$ds$a"]
      off4=m.topology.offset["$ds$a"]

      for k in 1:nd
        done=Int[];
        interim=Int[];
        kk=ind3[off3[k]:off3[k+1]-1]

        for h in ind1[off1[k]:off1[k+1]-1]
          for j in ind2[off2[h]:off2[h+1]-1]
            if in(j,done)
              continue;
            end

            jk=ind4[off4[j]:off4[j+1]-1]
            if subsetint(jk,kk)
              push!(interim,j);
              z+=1;
              t=true;
              push!(done,j);
            end
          end
        end

        sort!(interim);
        append!(i,interim);

        if t
          push!(o,z);
        end
        t=false;
      end
    end

    return i,o
  else
    error("Zum Berechnen der Inzidens d->d' durch meshIntersection muss gelten d>=d'!");
  end
end

#Funktion zum rekursiven Berechnen der Inzidenz d->d' durch Verwendung
#von meshTranspose und meshIntersection
function meshConnectivity!(m::mesh, d::Int, ds::Int)
  #Überprüfen, ob d->d´bereits vorhanden ist
  if !haskey(m.topology.incidence,"$d$ds")
    inc=Int[];
    off=Int[];

    #Durch Vergleich der Dimensionen miteinander, Bestimmen durch welche
    #Kombination von meshTranspose, meshIntersection und rekursiven Aufrufen von
    #meshConnectivity d->d' berechnet wird
    if d<ds
      meshConnectivity!(m,ds,d);
      inc,off=meshTranspose(m,d,ds);
    else
      if d==0 && ds==0
        dss=m.topology.D;
      else
        dss=0;
      end

      meshConnectivity!(m,d,dss);
      meshConnectivity!(m,dss,ds);
      inc,off=meshIntersection(m,d,ds,dss);
    end
    #Speichern der Inzidenz und des Offsets von d->d' in der Topologie des meshes m
    m.topology.incidence["$d$ds"]=inc;
    m.topology.offset["$d$ds"]=off;
  end
end

function completeMesh!(m)
  #jede mögliche Verbindung/Richtung wird erzeugt
    for i in 0:(m.topology.D)
        for j in 0:(m.topology.D)
          meshConnectivity!(m,i,j);
        end
    end
end

function checkSurr(m::mesh, d::Int, id::Int)
    #es müssen (fast) alle Verbindungen existieren, damit in jeder Dimension die Nachbarn
    #geprüft werden können
    completeMesh!(m)
    surr = Dict{String, Array{Int,1}}()

    for i in 0:(m.topology.D)
        #jede mögliche Dimension wird durchlaufen
        dir = "$d$i"
        #Feststellen des richtigen Bereichs im Inzidenzvektor
        rfirst = m.topology.offset[dir][id]
        rlast = m.topology.offset[dir][id+1]-1
        range = rfirst:rlast

        #Feststellen der Nachbarn über Inzidenzvektor
        current   = m.topology.incidence[dir][range]
        surr[dir] = current
    end

    return surr
end
