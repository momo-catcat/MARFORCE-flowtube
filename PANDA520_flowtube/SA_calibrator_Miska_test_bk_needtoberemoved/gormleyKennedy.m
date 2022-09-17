function n = gormleyKennedy( diffCoeff,tubeLength,lpm )

    e=pi*diffCoeff*tubeLength/(lpm/1000/60);
    
    n=zeros(size(e));
    
    for i=1:numel(e)
        if e(i)<0.02
            n(i)=1-2.56*e(i)^(2/3)+1.2*e(i)+0.177*e(i)^(4/3);
        else
            n(i)=0.819*exp(-3.657*e(i))+0.097*exp(-22.3*e(i))+0.032*exp(-57*e(i));
        end
    end

end

