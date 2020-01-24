macro solution(fieldname...)
           fields = [:($(Symbol(fieldname[i]))::Array{AbstractFloat,1}) for i=1:length(fieldname)]

           if isdefined(Main,:solution)
             if sort(collect(fieldname))!=sort(collect(fieldnames(solution)))
               error("Es existiert bereits ein struct solution. Um ein neues struct solution mit anderen field-Bezeichnungen zu erstellen, muss Julia mit exit() verlassen werden und das Programm neu ausgeführt werden.")
             else
               return nothing;
             end
           end

           s=esc(quote
               mutable struct solution
                   $(fields...)
                   solution()=new()
               end

               function createSolution(n::Int...)
                 s=solution();
                 k=fieldnames(solution);
                 length(k)!=length(n) && error("Die Anzahl der eingegebenen Längen müssen mit der Anzahl der Komponenten von solution übereinstimmen!")
                 for i in 1:length(n)
                   setfield!(s, k[i],zeros(n[i]));
                 end
                 return s;
               end

               function createSolution(v::Array{AbstractFloat,1}...)
                   s=solution();
                   k=fieldnames(solution);
                   length(k)!=length(v) && error("Die Anzahl der eingegebenen Matrizen müssen mit der Anzahl der Komponenten von solution übereinstimmen!")
                   for i in 1:length(v)
                     setfield!(s, k[i],v[i]);
                   end
                   return s;
               end

               import Base.+
               function +(s1::solution, s2::solution)
                 r=solution();
                 k=fieldnames(solution);
                 for i in 1:length(k)
                   a=getfield(s1, k[i])+getfield(s2, k[i]);
                   setfield!(r, k[i],a);
                 end
                 return r;
               end

               import Base.-
               function -(s1::solution, s2::solution)
                 r=solution();
                 k=fieldnames(solution);
                 for i in 1:length(k)
                   a=getfield(s1, k[i])-getfield(s2, k[i]);
                   setfield!(r, k[i],a);
                 end
                 return r;
               end

               import Base.*
               function *(val::AbstractFloat, s::solution)
                 r=solution();
                 k=fieldnames(solution);
                 for i in 1:length(k)
                   a=val*getfield(s, k[i]);
                   setfield!(r, k[i],a);
                 end
                 return r;
               end
           end)

           return s;
end
