function getReferenceBoundary(type::Symbol, mt::Int)
    if mt==4
        if type==:DG0

            cm=[0.5 0.0 0.0; 1.0 0.5 0.0; 0.5 1.0 0.0; 0.0 0.5 0.0];

        elseif type==:P1

            cm=[0.5 0.0 1.0 1.0 0.0 0.0; 1.0 0.5 0.0 1.0 1.0 0.0; 0.5 1.0 0.0 0.0 1.0 1.0; 0.0 0.5 1.0 0.0 0.0 1.0];

        #=
        elseif type==:P1x

            #cm=[0.5 0.0 0.0 0.0; 1.0 0.5 0.0 1.0; 0.5 1.0 0.0 0.0; 0.0 0.5 1.0 0.0];

        elseif type==:P1y

            cm=[0.5 0.0 1.0 0.0; 1.0 0.5 0.0 0.0; 0.5 1.0 0.0 1.0; 0.0 0.5 0.0 0.0];

        =#
        elseif type==:DG1

            cm=[0.5 0.0 0.0 0.0 0.0 0.0; 1.0 0.5 0.0 0.0 0.0 0.0; 0.5 1.0 0.0 0.0 0.0 0.0; 0.0 0.5 0.0 0.0 0.0 0.0];

        #=
        elseif type==:DG1x

            cm=[0.5 0.0 0.0 0.0; 1.0 0.5 0.0 0.0; 0.5 1.0 0.0 0.0; 0.0 0.5 0.0 0.0];

        elseif type==:DG1y

            cm=[0.5 0.0 0.0 0.0; 1.0 0.5 0.0 0.0; 0.5 1.0 0.0 0.0; 0.0 0.5 0.0 0.0];
        =#
        elseif type==:RT0

            cm=[0.5 0.0 1.0 0.0 0.0 0.0; 1.0 0.5 0.0 1.0 0.0 0.0; 0.5 1.0 0.0 0.0 1.0 0.0; 0.0 0.5 0.0 0.0 0.0 1.0];

        elseif type==:RT0B #Broken RT0

            cm=[0.5 0.0 0.0 0.0 0.0 0.0; 1.0 0.5 0.0 0.0 0.0 0.0; 0.5 1.0 0.0 0.0 0.0 0.0; 0.0 0.5 0.0 0.0 0.0 0.0];

        elseif type==:VecP1

            cm=[0.5 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0;
                1.0 0.5 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0;
                0.5 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0;
                0.0 0.5 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0];

        elseif type==:VecDG1

            cm=[0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                1.0 0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                0.5 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                0.0 0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];

        else
            error("Unzulässiger Finite-Elemente-Raum");
        end
    elseif mt==3
        if type==:DG0

            cm=[0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.0];

        elseif type==:RT0

            cm=[0.5 0.0 0.0 1.0 0.0;
                0.5 0.5 0.0 0.0 1.0;
                0.0 0.5 1.0 0.0 0.0];

        elseif type==:P1

            cm=[0.5 0.0 1.0 1.0 0.0;
                0.5 0.5 0.0 1.0 1.0;
                0.0 0.5 1.0 0.0 1.0];

        end
    end

    return cm
end
