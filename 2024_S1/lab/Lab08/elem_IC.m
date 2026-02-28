function fICEl = elem_IC(hl,Dl,Dl1,U0l,U0l1)
    slope = (U0l1 - U0l)/hl;
    fICEl = (slope*(Dl+Dl1)/2)*[1; -1];
end

