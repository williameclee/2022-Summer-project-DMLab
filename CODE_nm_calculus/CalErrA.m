function Err = CalErrA(F,G)
    Err = (mean(F,"all")-mean(G,"all"))^2+(std(F,[],"all")-std(G,[],"all"))^2;
end