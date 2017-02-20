
function q=qgsls(te,ts,L,w)
         
          % Q fuer den L-fachen standard linear solid:
            sumzQ=0;sumnQ=0;
            for l=1:L,
               d=1.0+w.*w*ts(l)*ts(l);
               sumzQ=((1.0+w.*w*te(l)*ts(l))./d)+sumzQ;
               sumnQ=(w*(te(l)-ts(l))./d)+sumnQ;
            end
          % Qualitaetsfaktor als Funktion der Frequenz fuer dem L-fachen SLS
           q=(1.0-L+sumzQ)./sumnQ;

