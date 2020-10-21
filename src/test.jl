using Combinatorics #for "permutations"

KommPot = [
    0 2 1 0 0 0 0 0 1 2 2;
    2 0 2 1 0 0 0 0 2 1 1;
    1 2 0 2 1 0 0 1 1 0 0;
    0 1 2 0 2 1 0 1 0 0 0;
    0 0 1 2 0 2 1 1 0 0 0;
    0 0 0 1 2 0 2 2 0 0 0;
    0 0 0 0 1 2 0 2 0 0 0;
    0 0 1 1 1 2 2 0 0 0 0;
    1 2 1 0 0 0 0 0 0 2 1;
    2 1 0 0 0 0 0 0 2 0 2;
    2 1 0 0 0 0 0 0 1 2 0
]

PersRel = [
    0 2 0 0 0 0 0 0 0 1 0	   ;
    2 0 0 -1 0 0 0 1 0 0 0    ;
    0 0 0 0 0 0 0 0 0 0 -2    ;
    0 -1 0 0 0 0 0 0 -2 0 -2  ;
    0 0 0 0 0 2 0 0 0 0 -2    ;
    0 0 0 0 2 0 -2 0 0 0 1    ;
    0 0 0 0 0 -2 0 0 1 0 -2   ;
    0 1 0 0 0 0 0 0 -2 0 0    ;
    0 0 0 -2 0 0 1 -2 0 2 0   ;
    1 0 0 0 0 0 0 0 2 0 0     ;
    0 0 -2 -2 -2 1 -2 0 0 0 0
]

# "All possible permutations" - https://en.wikibooks.org/wiki/Julia_for_MATLAB_Users/Core_Language/Mathematics#perms_All_possible_permutations
perms(a) = reverse(collect(permutations(a)))

function MaalFkt(Pl,Tr,Dst)
    N = length(Pl);
    MV = 0;
    for i = 1:N-1
     for j = i+1:N
     MV = MV + Tr[Pl[i],Pl[j]]*Dst[i,j];
     end
    end

    return MV
end

function OptPlacMaxMV(Tr,Dst)
    N = size(Tr)[1];
    m = factorial(N)
    a = collect(1:N)
    OptPl = collect(1:N)
    OptMV = -Inf; # Initialværdi - leder efterflg. efter løsn. med større værdi
    for i = 1:m
     Pl = nthperm(a,i);
     MV = MaalFkt(Pl,Tr,Dst);
     if MV > OptMV
        OptPl = Pl;
        OptMV = MV;
     end
    end
    return OptPl,OptMV
end

@time OptBordplan,OptScore = OptPlacMaxMV(PersRel,KommPot)