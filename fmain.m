clear all;
clc;

No = 0.01; %noise power  
do = 1; %reference distance       
a = 3;  %path loss component

global sou;
global rel;
global des;
global desr;

global d_sr;
global d_rd;
global d_sd;
global d_rrd;

global rel1;
global rel2;
global rel3;

k=1;
for i=1:6
    
  [sou, rel1, rel2, rel3]= with_relay_destination(i);

    for j=1:3

        %sou = [0 0];  %source
        %rel = [2 2];  %relay
        des = [8 8];  %destination of source
        desr = [10 10]; %destination of relay

        %fhandle=@optim_x_using_loop;
        if j==1
            rel=rel1;
        elseif j==2
            rel=rel2;
        elseif j==3
            rel=rel3;
        end
        s_r =  [sou;rel];   %source to relay
        r_d =  [rel;des];   %relay to source-destination
        s_d =  [sou;des];   %source to source-destination
        r_dr = [rel;desr]; %relay to relay-destination

        d_sr = pdist(s_r);   %source to relay
        d_rd = pdist(r_d);   %relay to source-destination
        d_sd = pdist(s_d);    %source to source-destination 
        d_rrd = pdist(r_dr);  %relay to relay-destination


        [x fval] = fminbnd(@optim_x_using_loop,0.0001,0.999)

        y(k,:) = [x -fval];
        uti(j,i) = -fval;

        % calculating data rate at optimum
         bw = 10, ps =1; pr=1;

         snr_sd1 = (ps/(No*bw)) * ( do / d_sd)^a;   %snr of source to destination channel
         r_sd1 = log2(1 + snr_sd1);             %data rate

            snr_sd  = (ps/(No*bw*(1-(x)))) * ( do / d_sd)^a;   
            snr_sr  = (ps/(No*bw*(1-(x)))) * ( do / d_sr)^a;    
            snr_rd  = (pr/(No*bw*(1-(x)))) * ( do / d_rd)^a;    
            snr_rrd  = (pr/(No*bw*(x)))* (do/d_rrd)^a;    
        R_srd(j,i) = (0.5* bw*(1-(x))) * min(log2(1+ snr_sr ), log2(1 + snr_sd  + snr_rd )); 
        r_rrd(j,i) = (bw*(x)) * log2(1 + snr_rrd ); %data rate from relay to relay-destination
        k=k+1;
    
    end

end


snr_sd1 = (ps/(No)) * ( do / d_sd)^a;
 r_sd1 = log2(1 + snr_sd1);
dataset_source={'No','source1','source2','source3','source4','source5','source6';...
    'Relay 1' R_srd(1) R_srd(4) R_srd(7) R_srd(10) R_srd(13) R_srd(16);...
    'Relay 2' R_srd(2) R_srd(5) R_srd(8) R_srd(11) R_srd(14) R_srd(17);...
    'Relay 3' R_srd(3) R_srd(6) R_srd(9) R_srd(12) R_srd(15) R_srd(18)};

dataset_relay={'No','source1','source2','source3','source4','source5','source6';...
    'Relay 1' r_rrd(1) r_rrd(4) r_rrd(7) r_rrd(10) r_rrd(13) r_rrd(16);...
    'Relay 2' r_rrd(2) r_rrd(5) r_rrd(8) r_rrd(11) r_rrd(14) r_rrd(17);...
    'Relay 3' r_rrd(3) r_rrd(6) r_rrd(9) r_rrd(12) r_rrd(15) r_rrd(18)};

 excel={'No','source1','source2','source3','source4','source5','source6';...
    'Relay 1' y(1) y(4) y(7) y(10) y(13) y(16);...
    'Relay 2' y(2) y(5) y(8) y(11) y(14) y(17);...
    'Relay 3' y(3) y(6) y(9) y(12) y(15) y(18)};

xlswrite('output.xlsx', excel) %  write data to untitle excell file
xlswrite('output1.xlsx', dataset_source)
xlswrite('output2.xlsx', dataset_relay)

R_srd;
r_rrd;
y;
uti;

n=6;
m=3;

S=[1 2 3 4 5 6];
R=[1 2 3];

for k=1:3
    

    [mS I] = max(R_srd, [], 1);
    [mR J] = max(r_rrd, [], 2);

    mSz=0.*R_srd;
    mRz=0.*r_rrd;

    for i=1:n

        mSz(I(i),i)=1;

    end

    for j=1:m

        mRz(j,J(j))=1;

    end

    mSRz=mSz.*mRz

    for i=1:m

        for j=1:n

            if mSRz(i,j)==1
                
                   
               %display(i+p)
                Source=S(j)
                Relay=R(i)

                R_srd = R_srd( [1:i-1,i+1:end] , : );
                R_srd = R_srd( : ,[1:j-1,j+1:end]  )
                r_rrd = r_rrd( [1:i-1,i+1:end] , : );
                r_rrd = r_rrd( : ,[1:j-1,j+1:end]  )
                
                S=S([1:j-1,j+1:end])
                R=R([1:i-1,i+1:end])
            end
            

        end

    end
    n=n-1;
    m=m-1;
end
        


%[mU K] = max(uti, [], 1)

   
   
        