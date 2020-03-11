function embed(comp::Symbol,degF::degF{1},cval::Array{Float64,1},compRec::Symbol,degFRec::degF{1},n::Int64)
    cEmbed=zeros(degFRec.numB);

    globalNum=Array{Int64,1}(undef,length(degF.phi));
    globalNumRec=Array{Int64,1}(undef,length(degFRec.phi));
    #println(" S: comp ",comp," compRec ",compRec)
    if (comp==:DG0 && compRec==:DG1) || (comp==:DG1 && compRec==:DG2)
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNumRec)
                for k in 1:length(globalNum)
                    cEmbed[globalNumRec[j]]+=cval[globalNum[k]];
                end
            end
        end
    elseif ((comp==:P1 || comp==:DG1) && compRec==:DG1) || ((comp==:P2 || comp==:DG2) && compRec==:DG2)
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNumRec)
                cEmbed[globalNumRec[j]]+=cval[globalNum[j]];
            end
        end
    elseif (comp==:DG1 && compRec==:P1) || (comp==:DG2 && compRec==:P2)
        z=zeros(degFRec.numB);
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
        z=zeros(degFRec.numB);
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

function embed(comp::Symbol,degF::degF{2},cval::Array{Float64,1},compRec::Symbol,degFRec::degF{2},n::Int64)
    cEmbed=zeros(degFRec.numB);

    globalNum=Array{Int64,1}(undef,size(degF.phi,2));
    globalNumRec=Array{Int64,1}(undef,size(degFRec.phi,2));

    #println(" V: comp ",comp," compRec ",compRec)
    if comp==:VecP1 && compRec==:VecDG1
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNum)
                cEmbed[globalNumRec[j]]+=cval[globalNum[j]];
            end
        end
    elseif comp==:VecP2 && compRec==:VecDG2
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNum)
                cEmbed[globalNumRec[j]]+=cval[globalNum[j]];
            end
        end
    elseif (comp==:RT0 || comp==:RT0B) && compRec==:VecDG1
        #h=[1,1,2,2,3,3,4,4];
        if length(globalNumRec)==8
            assign = [[[4],[1],[2],[1],[2],[3],[4],[3]]]
            factors = [[[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0]]]
        elseif length(globalNumRec)==6
            assign = [[[3],[1],[2,1],[1],[3],[2,3]],
                      [[3],[1,3],[3],[2],[1,2],[2]]]
            factors = [[[-1.0],[-1.0],[1.0,1.0],[-1.0],[-1.0],[1.0,1.0]],
                       [[-1.0],[1.0,1.0],[-1.0],[-1.0],[1.0,1.0],[-1.0]]]
        end

        cellshape=1;
        ncellshapes=length(assign);
        for i in 1:n
            cellshape=1+mod(i-1,ncellshapes)
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNumRec)
                for k in 1:length(assign[cellshape][j])
                    cEmbed[globalNumRec[j]]+=factors[cellshape][j][k]*cval[globalNum[assign[cellshape][j][k]]];
                end
            end
        end
        println("cval")
        println(cval)
        println("cembed")
        println(cEmbed)

#=
        #h=[1,1,2,2,3,3,4,4];
        if length(globalNumRec)==8
            h=[4,1,2,1,2,3,4,3]
        elseif length(globalNumRec)==6
            h=[3,1,2,1,3,2]
        end
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            for j in 1:length(globalNumRec)
                cEmbed[globalNumRec[j]]+=cval[globalNum[h[j]]];
            end
        end
=#
    elseif (comp==:RT1 || comp==:RT1B) && compRec==:VecDG2
        for i in 1:n
            l2g!(globalNum,degF,i);
            l2g!(globalNumRec,degFRec,i);
            cEmbed[globalNumRec[ 1]]=0.5*(cval[globalNum[1]]+cval[globalNum[3]])
            cEmbed[globalNumRec[ 2]]=0.5*(cval[globalNum[2]]+cval[globalNum[4]])
            cEmbed[globalNumRec[ 3]]=cval[globalNum[1]]
            cEmbed[globalNumRec[ 4]]=0.5*(cval[globalNum[5]]+cval[globalNum[6]])
            cEmbed[globalNumRec[ 5]]=0.5*(cval[globalNum[7]]+cval[globalNum[8]])
            cEmbed[globalNumRec[ 6]]=cval[globalNum[2]]
            cEmbed[globalNumRec[ 7]]=cval[globalNum[3]]
            cEmbed[globalNumRec[ 8]]=0.5*(cval[globalNum[9]]+cval[globalNum[10]])
            cEmbed[globalNumRec[ 9]]=0.5*(cval[globalNum[11]]+cval[globalNum[2]])
            cEmbed[globalNumRec[10]]=cval[globalNum[4]]
            cEmbed[globalNumRec[11]]=cval[globalNum[11]]
            cEmbed[globalNumRec[12]]=cval[globalNum[5]]
            cEmbed[globalNumRec[13]]=cval[globalNum[7]]
            cEmbed[globalNumRec[14]]=cval[globalNum[6]]
            cEmbed[globalNumRec[15]]=cval[globalNum[8]]
            cEmbed[globalNumRec[16]]=cval[globalNum[10]]
            cEmbed[globalNumRec[17]]=cval[globalNum[12]]
            cEmbed[globalNumRec[18]]=cval[globalNum[9]]
        end
    else
        error("Entsprechende embed-Funktion fehlt.")
    end

    return cEmbed;
end
