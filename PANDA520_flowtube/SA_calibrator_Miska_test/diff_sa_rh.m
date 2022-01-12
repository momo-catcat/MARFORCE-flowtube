function D = diff_sa_rh(T,rh)
%DIFF_SA_RH Diffusion coefficient of sulfuric acid that takes rh into
%account. T in Kelvins (can be a vector), rh between 0 and 3 (not percents).


rh = min(rh,3);

D = 1.8e-9*T.^1.5/(1+.2876*rh^.5643); %Fitted to measured data of Hanson&Eisele (2000)

end

