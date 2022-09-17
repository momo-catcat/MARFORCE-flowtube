function [ tulos ] = vappresw( T )
%veden saturaatiohöyrynpaine (Pa) lämpötilassa T (K)
%Preining 1981

tulos = exp(77.34491296-7235.424651./T-8.2.*log(T)+0.0057113.*T);

end

