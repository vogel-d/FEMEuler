macro diagnostic(fieldname...)
           fields = [:($(Symbol(fieldname[i]))::Array{Float64,1}) for i=1:length(fieldname)]

           if isdefined(Main,:diagnostic)
             if sort(collect(fieldname))!=sort(collect(fieldnames(diagnostic)))
               error("Es existiert bereits ein struct diagnostic. Um ein neues struct diagnostic mit anderen field-Bezeichnungen zu erstellen, muss Julia mit exit() verlassen werden und das Programm neu ausgeführt werden.")
             else
               return nothing;
             end
           end

           s=esc(quote
               mutable struct diagnostic
                   $(fields...)
                   diagnostic()=new()
               end

               function createDiagnostic(n::Int64...)
                 s=diagnostic();
                 k=fieldnames(diagnostic);
                 length(k)!=length(n) && error("Die Anzahl der eingegebenen Längen müssen mit der Anzahl der Komponenten von diagnostic übereinstimmen!")
                 for i in 1:length(n)
                   setfield!(s, k[i],zeros(n[i]));
                 end
                 return s;
               end

               function createDiagnostic(v::Array{Float64,1}...)
                   s=diagnostic();
                   k=fieldnames(diagnostic);
                   length(k)!=length(v) && error("Die Anzahl der eingegebenen Matrizen müssen mit der Anzahl der Komponenten von diagnostic übereinstimmen!")
                   for i in 1:length(v)
                     setfield!(s, k[i],v[i]);
                   end
                   return s;
               end

               import Base.+
               function +(s1::diagnostic, s2::diagnostic)
                 r=diagnostic();
                 k=fieldnames(diagnostic);
                 for i in 1:length(k)
                   a=getfield(s1, k[i])+getfield(s2, k[i]);
                   setfield!(r, k[i],a);
                 end
                 return r;
               end

               import Base.-
               function -(s1::diagnostic, s2::diagnostic)
                 r=diagnostic();
                 k=fieldnames(diagnostic);
                 for i in 1:length(k)
                   a=getfield(s1, k[i])-getfield(s2, k[i]);
                   setfield!(r, k[i],a);
                 end
                 return r;
               end

               import Base.*
               function *(val::Float64, s::diagnostic)
                 r=diagnostic();
                 k=fieldnames(diagnostic);
                 for i in 1:length(k)
                   a=val*getfield(s, k[i]);
                   setfield!(r, k[i],a);
                 end
                 return r;
               end
           end)

           return s;
end
