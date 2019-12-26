using SymPy;

function getCoefficients(p::Function)
    @vars x y;
    f=expand(p(x,y));
    coeff=f.as_coefficients_dict();
    n=real(maximum([degree(f,x),degree(f,y)]));
    res=zeros(Float64,n+1,n+1)
    for k in keys(coeff)
        nx=real(degree(k,x));
        ny=real(degree(k,y));
        res[ny+1,nx+1]=float(coeff[k]);
    end
    print("[")
    for i in 1:size(res,2)
        for j in 1:size(res,1)
            if i==1 && j==1
                print("$(res[i,j])")
            else
                print(" $(res[i,j])")
            end
        end
        i!=size(res,2) && print(";")
    end
    print("]")
    return res;
end
