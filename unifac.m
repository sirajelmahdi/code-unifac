%variable definition et calcul des parametre 
nb=input('entrer le nombre d elemnet ');
temp=input('entrer la temperature');
coeffexp=input('entrer la valeur de coeff diffusion experimental ');
sommteta=0;
sommephi=0;
for i=1:(nb^2)
   if i<=nb 
       if i==1
x(i,1)=input('entrer la fraction solute en 1er lieu a  ');

r(i,1)=input('entrer la r 1 de ssolute a  ');

q(i,1)=input('entrer la q 1 de solute  a  ');

D0(i,1)=input('diffusion 1 a faible concentration dab ');
A(i,1)=input('les valeure de a11 solute'); 
to(i,1)=exp((-A(i,1))/temp);
       end
       if i==2
      x(i,1)=input('entrer la fraction solvant b ');

    r(i,1)=input('entrer la r  solvant b ');

    q(i,1)=input('entrer la q  sollute b  ');

    D0(i,1)=input('diffusion 2 a faible concentration dba ');
    A(i,1)=input('les valeure de a12 ab'); 
    to(i,1)=exp((-A(i,1))/temp);
       end
       if i<=2
lamda(i,1)=(r(i,1))^(1/3);
sommteta=sommteta+(x(i,1)*(q(i,1)));
sommephi =sommephi+(x(i,1)*(lamda(i,1)));
       end
   end
   
   if i==3
 A(i,1)=input('les valeure de a21 ba'); 
 to(i,1)=exp((-A(i,1))/temp);
   end 
   if i==4
      A(i,1)=input('les valeure de a22 bb'); 
      to(i,1)=exp((-A(i,1))/temp);
   end


end
for i=1:(nb^2)
if i<=nb
teta(i,1)=(x(i,1)*(q(i,1)))/sommteta;
phiie(i,1)=(x(i,1)*(lamda(i,1)))/sommephi;
end


if i==3
somab=to(i-1,1)*teta(i-2,1)+to(i+1,1)*teta(i-1,1)
tetaab=to(i-1,1)*teta(i-2,1)/somab;
tetabb=to(i+1,1)*teta(i-1,1)/somab;
end

if i==4
somba=to(i-3,1)*teta(i-3,1)+to(i-1,1)*teta(i-2,1)
tetaba=to(i-1,1)*teta(i-2,1)/somba;
tetaaa=(to(i-3,1)*teta(i-3,1))/somba
end

end
 text   =     'le 11 cest le solute a '; 
 text1  =     'le 12 cest ab ';
 text2  =     'le 21 cest ba ';
 text3  =     'le 22 cest le solvant b';    
disp(text);disp(text1);disp(text2);disp(text3);
%partie de calcule du coeffcient de 1 er itteration ou A12 =A21=1
part1=(x(2,1)*log(D0(1,1)))+(x(1,1)*(D0(2,1)))+2*((x(1,1)*log((x(1,1)/phiie(1,1))))+x(2,1)*log(x(2,1)/phiie(2,1)));
part2=2*x(2,1)*x(1,1)*(((phiie(1,1)/x(1,1))*(1-(lamda(1,1)/lamda(2,1))))+((phiie(2,1)/x(2,1))*(1-(lamda(2,1)/lamda(1,1)))));
part3=x(2,1)*q(1,1)*(((1-tetaba^2)*log(to(3,1))+(1-tetabb^2)*log(to(2,1))));
part4=x(1,1)*q(2,1)*(((1-tetaab^2)*log(to(2,1))+(1-tetaaa^2)*log(to(3,1))));
coefdiffth=exp(part1+part2+part3+part4)
k=1;;ind=0;inds=0;
A(1,1)=1;
A(2,1)=1;
A(3,1)=1;
A(4,1)=1;
A(5,1)=1;
A(6,1)=1;
A(7,1)=1;
A(8,1)=1;
%code specefique au cas du coeff exper inferieur du celle calculer en
%premiere itteration ===>pour avoir convergence
   if coefdiffth > coeffexp
        while coefdiffth-coeffexp < 0.0009 & delta == 1
            %code specif du variation de a12 et a21 en 4 variation
            A(1,1)=A(1,1)-0.1;%le cas a12 -0.1
            A(2,1)=A(2,1)-0.1;%le cas a21 -0.1
            A(3,1)=A(3,1)+0.1;%le cas a12 +0.1
            A(4,1)=A(4,1)+0.1;%le cas a21 +0.1
            A(5,1)=A(5,1)+0.1;%le cas a12 +0.1
            A(6,1)=A(6,1)-0.1;%le cas a21 -0.1
            A(7,1)=A(7,1)-0.1;%le cas a12 -0.1
            A(8,1)=A(8,1)+0.1;%le cas a21 +0.1
            k=k+1
            %code pour calcule le coeff theorique  pour chaque cas de
            %variation de a12 a21
                    for i=1:8
                            if mod(i,2)~= 0 
                            abto(i,1)=exp(-A(i,1)/temp);
                            bato(i,1)=exp(-A(i+1,1)/temp);
                            abteta(i,1)=teta(1,1)*abto(i,1)/(teta(1,1)*abto(i,1)+teta(2,1));
                            bateta(i,1)=teta(2,1)*bato(i,1)/(teta(2,1)*bato(i,1)+teta(1,1));
                            aateta(i,1)=teta(1,1)/(teta(1,1)+teta(2,1)*bato(i,1));
                            bbteta(i,1)=teta(2,1)/(teta(2,1)+teta(1,1)*abto(i,1));
                            part(1,1)=x(2,1)*q(1,1)*((1-bateta(i,1)^2)*log(bato(i,1))+(1-bbteta(i,1)^2)*abto(i,1)*log(abto(i,1)));
                            part(2,1)=x(1,1)*q(2,1)*((1-abteta(i,1)^2)*log(abto(i,1))+(1-aateta(i,1)^2)*bato(i,1)*log(bato(i,1)));
                            cas(i,1)=exp(part(1,1)+part(2,1)+part1+part2);
                            inc(k,i)=cas(i,1);
                            indab(k,i)=A(i,1);
                            indba(k,i)=A(i+1,1);
                            incx(k,i)=abto(i,1);
                            incy(k,i)=bato(i,1);
                                    if i==1
                                    coed(1,1)=cas(i,1);
                                    end
                                        if i==3
                                        coed(2,1)=cas(i,1);
                                        end
                                            if i==5
                                            coed(3,1)=cas(i,1);
                                            end
                                                if i==7
                                                coed(4,1)=cas(i,1);
                                                end
                            end

                    end
                    %detection  valeur plus proche a coef exp par
                    %determination du minimum des coef theo  pour le cas du
                    %coef theo 0 inf du coef exp
                            minimum(k,1)=min(coed);
%                             coefdiffth=minimum(k,1);
                                 %verificat que le minimum est entre coef
                                 %theo 0 et coefexp et le codition d arret
                                
                                     if minimum(k,1) > coeffexp
                                         coefdiffth=minimum(k,1);
                                         j=j+1
                                         delta = 1
                                     else
                                         delta=0
                                     end
                                 
            

        end
   end
   %code specefique si cas ou le coefficient experimentale est superieure
   %de celle calculer en premiere itteraion
   %le meme  algorithme et raisonnement pour le cas ou le coefficient
   %experimentale est superieure coef theo 0
    if coefdiffth < coeffexp
        while coefdiffth-coeffexp < 0.0009 & delta == 1
            A(1,1)=A(1,1)-0.1;
            A(2,1)=A(2,1)-0.1;
            A(3,1)=A(3,1)+0.1;
            A(4,1)=A(4,1)+0.1;
            A(5,1)=A(5,1)+0.1;
            A(6,1)=A(6,1)-0.1;
            A(7,1)=A(7,1)-0.1;
            A(8,1)=A(8,1)+0.1;
            k=k+1
                    for i=1:8
                            if mod(i,2)~= 0 
                            abto(i,1)=exp(-A(i,1)/temp);
                            bato(i,1)=exp(-A(i+1,1)/temp);
                            abteta(i,1)=teta(1,1)*abto(i,1)/(teta(1,1)*abto(i,1)+teta(2,1));
                            bateta(i,1)=teta(2,1)*bato(i,1)/(teta(2,1)*bato(i,1)+teta(1,1));
                            aateta(i,1)=teta(1,1)/(teta(1,1)+teta(2,1)*bato(i,1));
                            bbteta(i,1)=teta(2,1)/(teta(2,1)+teta(1,1)*abto(i,1));
                            part(1,1)=x(2,1)*q(1,1)*((1-bateta(i,1)^2)*log(bato(i,1))+(1-bbteta(i,1)^2)*abto(i,1)*log(abto(i,1)));
                            part(2,1)=x(1,1)*q(2,1)*((1-abteta(i,1)^2)*log(abto(i,1))+(1-aateta(i,1)^2)*bato(i,1)*log(bato(i,1)));
                            cas(i,1)=exp(part(1,1)+part(2,1)+part1+part2);
                            inc(k,i)=cas(i,1);
                            indab(k,i)=A(i,1);
                            indba(k,i)=A(i+1,1);
                            incx(k,i)=abto(i,1);
                            incy(k,i)=bato(i,1);
                                    if i==1
                                    coed(1,1)=cas(i,1);
                                    end
                                        if i==3
                                        coed(2,1)=cas(i,1);
                                        end
                                            if i==5
                                            coed(3,1)=cas(i,1);
                                            end
                                                if i==7
                                                coed(4,1)=cas(i,1);
                                                end
                            end

                    end
                            maximum(k,1)=max(coed);
%                             coefdiffth=minimum(k,1);
                                 
                                     if maximum(k,1) < coeffexp
                                         coefdiffth=maximum(k,1);
                                         j=j+1
                                         delta = 1
                                     else
                                         delta=0
                                     end
                                 
            

        end
    end
  %localisation du colone et ligne de solution 
  for i = 1:k
      for e=1:7
      if coefdiffth==inc(i,e)
          col=i
          lin=e
      end
      end
  end
  %determination du A12 et A21 apartir des tableau du localisation
  Aab=indab(col,lin)
  Aba=indba(col,lin)
  
  
   %verification 
        abto(6,1)=exp(-Aab/temp);
        bato(6,1)=exp(-Aba/temp);
        abteta(6,1)=teta(1,1)*abto(6,1)/(teta(1,1)*abto(6,1)+teta(2,1));
        bateta(6,1)=teta(2,1)*bato(6,1)/(teta(2,1)*bato(6,1)+teta(1,1));
        aateta(6,1)=teta(1,1)/(teta(1,1)+teta(2,1)*bato(6,1));
        bbteta(6,1)=teta(2,1)/(teta(2,1)+teta(1,1)*abto(6,1));
        part(1,1)=x(2,1)*q(1,1)*((1-bateta(6,1)^2)*log(bato(6,1))+(1-bbteta(6,1)^2)*abto(6,1)*log(abto(6,1)));
        part(2,1)=x(1,1)*q(2,1)*((1-abteta(6,1)^2)*log(abto(6,1))+(1-aateta(6,1)^2)*bato(6,1)*log(bato(6,1)));
        verf=exp(part(1,1)+part(2,1)+part1+part2)
 