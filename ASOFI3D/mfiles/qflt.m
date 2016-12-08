% This function computes the difference between a
% frequency indepent quality factor and Q as function of
% relaxation frequencies and tau (written for optimization with leastsq).

function delta=qstd(x)

         global L w Qf1 Qf2
         
         fl=x(1:L);
         t=x(L+1);

         % computing telaxation times and relaxation frequencies
           ts=1./(2*pi*fl);

          % Q for a generalized standard linear solid:
            sumzQ=0;sumnQ=0;
            for l=1:L,
               d=1+w.*w*ts(l)*ts(l);
               sumnQ=(w*ts(l)*t./d)+sumnQ;
               sumzQ=((w.*w*ts(l)*ts(l)*t)./d)+sumzQ;
            end

           Qf2=(1+sumzQ)./sumnQ;

          delta=(Qf2-Qf1);

