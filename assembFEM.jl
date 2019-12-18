#Funktion zum Berechnen der Masse- und Steifigkeitsmatrizen
#In cond gibt der erste Wert die RB für die linke &rechte Seite des Meshes an und der zweite Wert für die obere & untere Seite
function assembFEM!(p::femProblem)
    assembMass!(p);
    assembStiff!(p);
    return nothing;
end
