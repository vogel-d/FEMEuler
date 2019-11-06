function embed(comp::Symbol,degF::degF{3},cval::Array{Float64,1},compRec::Symbol,degFRec::degF{3},n::Int64)
    cEmbed=zeros(size(degFRec.coordinates,2));

    globalNum=Array{Int64,1}(undef,size(degF.phi,3));
    globalNumRec=Array{Int64,1}(undef,size(degFRec.phi,3));
    if comp==:DG0 && compRec==:DG1
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);

            for j in 1:length(globalNumRec)
                for k in 1:length(globalNum)
                    cEmbed[globalNumRec[j]]+=cval[globalNum[k]];
                end
            end
        end
    elseif ((comp==:P1 || comp==:DG1) && compRec==:DG1)
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNumRec)
                cEmbed[globalNumRec[j]]+=cval[globalNum[j]];
            end
        end
    elseif comp==:DG1 && compRec==:P1
        z=zeros(size(degFRec.coordinates,2));
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNumRec)
                cEmbed[globalNumRec[j]]+=cval[globalNum[j]];
                z[globalNumRec[j]]+=1;
            end
        end
        for i in 1:length(z)
            if z[i]==0.0
                z[i]=1.0
            end
        end
        cEmbed=cEmbed./z;
    elseif (comp==:P1y || comp==:DG1y )&& compRec==:P1
        z=zeros(size(degFRec.coordinates,2));
        h=[1,1,2,2];
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNumRec)
                cEmbed[globalNumRec[j]]+=cval[globalNum[h[j]]];
                z[globalNumRec[j]]+=1;
            end
        end
        for i in 1:length(z)
            if z[i]==0.0
                z[i]=1.0
            end
        end
        cEmbed=cEmbed./z;
    else
        error("Entsprechende embed-Funktion fehlt.")
    end

    return cEmbed;
end

function embed(comp::Symbol,degF::degF{4},cval::Array{Float64,1},compRec::Symbol,degFRec::degF{4},n::Int64)
    cEmbed=zeros(size(degFRec.coordinates,2));

    globalNum=Array{Int64,1}(undef,size(degF.phi,4));
    globalNumRec=Array{Int64,1}(undef,size(degFRec.phi,4));

    if comp==:VecP1 && compRec==:VecDG1
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:8
                cEmbed[globalNumRec[j]]+=cval[globalNum[j]];
            end
        end
    elseif (comp==:RT0 || comp==:RT0B) && compRec==:VecDG1
        h=[1,1,2,2,3,3,4,4];
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:8
                cEmbed[globalNumRec[j]]+=cval[globalNum[h[j]]];
            end
        end
    else
        error("Entsprechende embed-Funktion fehlt.")
    end

    return cEmbed;
end
