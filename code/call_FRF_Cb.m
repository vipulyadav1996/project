clc;
clear all;
wf = 0.1:0.4:30;
hold on
colour = ['r','g','b','y' ,'k'];
j = 1;p=1000;% p=10000
for F0 =[1000] %[0.42]

    for k2 = 1005:2:1015 %2.347 is optimum linear parameters.
                   %for secondary non-linearity k2opt=2.19,c2opt=0.34
    for c2= 25.5
        i = 1;
    for w=wf          %0.05:0.01:3.35
        [ssv(i),w_n]=FRF_Cb(w,F0,k2,c2,'n');
        ssv(i)=ssv(i);
        i = i+1;
    end
    pvalue=max(ssv)
    if pvalue < p
            p = pvalue;
            k2opt = k2
            c2opt = c2
    end
    %plot(wf./(w_n),ssv)
    j = j+1;
end
    end
end
legend('-DynamicLegend');