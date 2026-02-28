function KEl = elem_stiff(hl, Dl, Dl1)
    KEl = ((Dl+Dl1)/(2*hl))*[1, -1; -1, 1];
end

