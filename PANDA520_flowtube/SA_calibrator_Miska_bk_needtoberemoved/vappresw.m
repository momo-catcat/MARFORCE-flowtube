function [ tulos ] = vappresw( T )
%veden saturaatioh�yrynpaine (Pa) l�mp�tilassa T (K)
%Preining 1981

tulos = exp(77.34491296-7235.424651./T-8.2.*log(T)+0.0057113.*T);

end

