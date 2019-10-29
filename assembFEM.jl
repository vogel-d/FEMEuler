#Funktion zum Berechnen der Masse- und Steifigkeitsmatrizen
#In cond gibt der erste Wert die RB für die linke &rechte Seite des Meshes an und der zweite Wert für die obere & untere Seite
function assembFEM!(p::femProblem, cond::Tuple{Symbol,Symbol})
    generateEquals!(p,cond);
    if cond[1]==cond[2]
        if cond[1]==:periodic
            generatePeriodicBoundary!(p::femProblem);
        else
            p.degFBoundary=deepcopy(p.degF);
        end
    else
        generateMixedBoundary!(p::femProblem);
    end
    assembMass!(p);
    assembStiff!(p);
    return nothing;
end
