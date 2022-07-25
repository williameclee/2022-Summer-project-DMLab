function Err = CalErrF(F,G)
    rho = corrcoef(F,G);
    Err = 2*(1-rho(2))*std(F,[],"all")*std(G,[],"all");
end