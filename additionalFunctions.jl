import Base.+
function +(a::SparseVector{Float64,Int64}, b::Array{Float64,2})
    rows=collect(1:length(b));
    vals=Float64[a[i]+b[i] for i in 1:length(b)];
    return sparsevec(rows,vals)
end

#Funktion zum Testen, ob die int Elemente des Arrays k
#im Array g enthalten sind
function subsetint(k::Array{Int64,1},g::Array{Int64,1})
  t=true;
  for h in 1:length(k)
    if !any(isequal(k[h]), g)
      t=false;
    end
  end
  return t;
end

#Funktion zum vollständigen printen einer Matrix beliebiger Größe
#falls vollständiges printen nicht notwendig: display(Matrix)
function printMatrix(m::Array)
    maxPreDigits=floor(maximum(m));
    maxPreDigits=length(String("$maxPreDigits"));
    for i in 1:size(m,1)
        print(" "^(1+maxPreDigits-length(String("$(floor(m[i,1]))"))));
        for j in 1:size(m,2)
            current=round(m[i,j];digits=4);
            print(current);
            if j<size(m,2)
                next=m[i,j+1];
                next=floor(next);
                current=round(current-floor(current);digits=4);
                ws=9+maxPreDigits-length(String("$current"*"$next"));
                print(" "^ws);
            end
        end
        print("\n");
    end
end

function mean(v::Array{Float64, 1})
    return 1/length(v)*sum(v);
end
