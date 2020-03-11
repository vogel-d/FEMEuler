function getSpace(femType::Symbol)

    H1=[:DG0,:DG1,:DG2,:P1,:P2];
    H1div=[:RT0,:RT1,:RT0B,:RT1B];
    H1xH1=[:VecDG1,:VecP1];

    if in(femType,H1)
        return :H1
    elseif in(femType,H1div)
        return :H1div
    elseif in(femType,H1xH1)
        return :H1xH1
    else
        error("Bitte für $femType Raum spezifizieren.")
    end
end
