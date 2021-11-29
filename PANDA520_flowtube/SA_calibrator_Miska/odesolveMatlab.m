function cout = odesolveMatlab(timesteps,Zgrid,Rgrid,dt,kSO2pOH,kOHpHO2,kOHpOH,kSO3p2H2O,kHSO3pO2,O2conc,H2Oconc,SO2conc,D,R,L,Q,cc)
    % This is a Matlab version of the odesolve.c (University of Helsinki)
    % made and corrected by Miska Olin (Tampere University) 


    c=zeros(5,Rgrid,Zgrid);
    for i=1:5
        c(i,:,:)=cc{i};
    end

    initc=c;
    
    dz=L/(Zgrid-1);
    dr=R/(Rgrid-1);
    
    for m=1:timesteps
        for l=1:5
            for i=2:Zgrid-1
                for j=2:Rgrid/2
                    
                    r = (Rgrid - j) * dr;
                    r = abs(2. * r - R);
                    if i == 2
                    end
%                     disp(l),disp(i),disp(j),disp(r), disp(R)
                    term1 = D(l) * (1. / r * (initc(l,j,i) - initc(l,j - 1,i)) / (-2*dr) + (initc(l,j + 1,i) - 2. * initc(l,j,i) + initc(l,j - 1,i)) / (4 * dr * dr) + (initc(l,j,i + 1) - 2. * initc(l,j,i) + initc(l,j,i - 1)) / (dz * dz));
                    term2 = (2. * Q) / (pi * R^2) * (1. - r^2 / (R^2)) * (initc(l,j,i) - initc(l,j,i - 1)) / dz;
                    
                    
                    if l==1
                        term3=kSO2pOH * SO2conc * initc(5,j,i) - kHSO3pO2 * initc(1,j,i) * O2conc;
                    elseif l==2
                        term3=kHSO3pO2 * O2conc * initc(1,j,i) - kSO3p2H2O * H2Oconc * H2Oconc * initc(2,j,i);
                    elseif l==3
                        term3=kHSO3pO2 * O2conc * initc(1,j,i) - kOHpHO2 * initc(3,j,i) * initc(5,j,i);
                    elseif l==4
                        term3=kSO3p2H2O * H2Oconc * H2Oconc * initc(2,j,i);
                    elseif l==5
                        term3=-kSO2pOH * SO2conc * initc(5,j,i) - kOHpHO2 * initc(3,j,i) * initc(5,j,i) - 2*kOHpOH * initc(5,j,i) * initc(5,j,i);
                    end
                    
                    c(l,j,i) = dt * (term1 - term2 + term3) + initc(l,j,i);
                end
                
                c(l,Rgrid / 2 + 1,i) = c(l,Rgrid / 2,i);
            end
            
            i=Zgrid;
            for j=2:Rgrid/2

                r = (Rgrid - j) * dr;
                r = abs(2. * r - R);
                term1 = D(l) * (1. / r * (initc(l,j,i) - initc(l,j - 1,i)) / (-2*dr) + (initc(l,j + 1,i) - 2. * initc(l,j,i) + initc(l,j - 1,i)) / (4* dr * dr) + (initc(l,j,i) - 2. * initc(l,j,i-1) + initc(l,j,i - 2)) / (dz * dz));
                term2 = (2. * Q) / (pi * R^2) * (1. - r^2 / (R^2)) * (initc(l,j,i) - initc(l,j,i - 1)) / dz;

                if l==1
                    term3=kSO2pOH * SO2conc * initc(5,j,i) - kHSO3pO2 * initc(1,j,i) * O2conc;
                elseif l==2
                    term3=kHSO3pO2 * O2conc * initc(1,j,i) - kSO3p2H2O * H2Oconc * H2Oconc * initc(2,j,i);
                elseif l==3
                    term3=kHSO3pO2 * O2conc * initc(1,j,i) - kOHpHO2 * initc(3,j,i) * initc(5,j,i);
                elseif l==4
                    term3=kSO3p2H2O * H2Oconc * H2Oconc * initc(2,j,i);
                elseif l==5
                    term3=-kSO2pOH * SO2conc * initc(5,j,i) - kOHpHO2 * initc(3,j,i) * initc(5,j,i) - 2*kOHpOH * initc(5,j,i) * initc(5,j,i);
                end

                c(l,j,i) = dt * (term1 - term2 + term3) + initc(l,j,i);

            end

            c(l,Rgrid / 2 + 1,i) = c(l,Rgrid / 2,i);
        end
        
        initc=c;
    end
            
    
    cout = cell(1, 5);
    cout = cellfun(@(x)zeros(Rgrid, Zgrid), cout, 'UniformOutput', false);
    for i=1:5
        cout{i}=squeeze(c(i,:,:));
    end
end

