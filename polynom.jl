function polynom(a::Array{T},x::T,y::T) where T
    sol=zero(T)
    for j in 1:size(a,2)
        for i in 1:size(a,1)
            sol+=a[i,j]*x^(j-1)*y^(i-1);
        end
    end
    return sol;
end

#=
Schema:
        1     x      x^2       x^3       x^4  ...

1       1     x      x^2       x^3       x^4     ...
y       y     xy     x^2y      x^3y      x^4y    ...
y^2     y^2   xy^2   x^2y^2    x^3y^2    x^4y^2  ...
y^3     y^3   xy^3   x^2y^3    x^3y^3    x^4y^3  ...
y^4     y^4   xy^4   x^2y^4    x^3y^4    x^4y^4  ...
.
.
.
=#
