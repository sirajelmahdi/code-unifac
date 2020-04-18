clear all clc 
%% definition des variables par utilisateur  et calcule des parametres 
% nb=input('entrer le nombre d elemnet ');
% temp=input('entrer la temperature');
% coeffexp=input('entrer la valeur de coeff diffusion experimental ');
% sommteta=0;
% sommephi=0;
% for i=1:(nb^2)
%    if i<=nb 
%        if i==1
% x(i,1)=input('entrer la fraction solute en 1er lieu a  ');
% 
% r(i,1)=input('entrer la r 1 de ssolute a  ');
% 
% q(i,1)=input('entrer la q 1 de solute  a  ');
% 
% D0(i,1)=input('diffusion 1 a faible concentration dab '); A(i,1)=input('les valeure de a11 aa'); 
%
%
%        end
%        if i==2
%       x(i,1)=input('entrer la fraction solvant b ');
% 
%     r(i,1)=input('entrer la r  solvant b ');
% 
%     q(i,1)=input('entrer la q  sollute b  ');
% 
%     D0(i,1)=input('diffusion 2 a faible concentration dba ');A(i,1)=input('les valeure de a12 ab')
%    
%        end
%        if i<=2
% lamda(i,1)=(r(i,1))^(1/3);
% 
%        end
%    end
%    
%    if i==3
%  A(i,1)=input('les valeure de a21 ba'); 
%
%    end 
%    if i==4
%       A(i,1)=input('les valeure de a22 bb'); 
%    
%    end
% 
% 
% end
% %
% % end
%  text   =     'le 11 cest le solute a '; 
%  text1  =     'le 12 cest ab ';
%  text2  =     'le 21 cest ba ';
%  text3  =     'le 22 cest le solvant b'; 
%  a12=A(2,1);
%  a21=A(3,1);
% disp(text);disp(text1);disp(text2);disp(text3);
% %partie de calcule du coeffcient de 1 er itteration ou A12 =A21=1
%definition des valeur des variables 
x=[0.35;0.65];D0=[0.0000210;0.0000267];
A=[0;928.051753999 ;928.051753999;0];err=0.000000000001 ;% err c est le pas de variation  
sa=A(2,1);sb=A(3,1);
temp=313;coeffexp=1.33*10^-5
q=[1.432;1.4];
r=[1.4311;0.92];jj=0;
%caclcul des parametre initiale 
lamda=[r(1,1)^(1/3);r(2,1)^(1/3)];
to=[exp(-A(1,1)/temp);exp(-A(2,1)/temp);exp(-A(3,1)/temp);exp(-A(4,1)/temp)];
%a
teta(1,1)=(x(1,1)*(q(1,1)))/(((x(1,1)*(q(1,1)))+(x(2,1)*(q(2,1)))));
phiie(1,1)=(x(1,1)*(lamda(1,1)))/((x(1,1)*(lamda(1,1)))+(x(2,1)*(lamda(2,1))));

%b
teta(2,1)=(x(2,1)*(q(2,1)))/(((x(1,1)*(q(1,1)))+(x(2,1)*(q(2,1)))));
phiie(2,1)=(x(2,1)*(lamda(2,1)))/((x(1,1)*(lamda(1,1)))+(x(2,1)*(lamda(2,1))));
tetabb=teta(2,1)/(to(2,1)*teta(1,1)+teta(2,1));
%ab
tetaab=to(2,1)*teta(1,1)/(to(2,1)*teta(1,1)+teta(2,1));
tetaaa=teta(1,1)/((to(3,1)*teta(2,1))+teta(1,1));
%ba
tetaba=to(3,1)*teta(2,1)/(to(3,1)*teta(2,1)+teta(1,1));
tetabb=teta(2,1)/(to(2,1)*teta(1,1)+teta(2,1));
%%section de calcul de 1 er coeff  en utilisant les valeur initialiser par l utilisateur 
part1=(x(2,1)*log(D0(1,1)))+(x(1,1)*(D0(2,1)))+2*((x(1,1)*log((x(1,1)/phiie(1,1))))+x(2,1)*log(x(2,1)/phiie(2,1)));
part2=2*x(2,1)*x(1,1)*(((phiie(1,1)/x(1,1))*(1-(lamda(1,1)/lamda(2,1))))+((phiie(2,1)/x(2,1))*(1-(lamda(2,1)/lamda(1,1)))));
part3=x(2,1)*q(1,1)*(((1-tetaba^2)*log(to(3,1))+(1-tetabb^2)*to(2,1)*log(to(2,1))));
part4=x(1,1)*q(2,1)*(((1-tetaab^2)*log(to(2,1))+(1-tetaaa^2)*to(3,1)*log(to(3,1))));
coef00=exp(part1+part2+part3+part4) % le coeff calculer pour a12 a21 initialiser par l'utilisateur
%% 	code spécifique à la résolution par fonction fminsearch
format longEng
uniff = @(a) abs(exp((x(2,1)*log(D0(1,1)))+(x(1,1)*(D0(2,1)))+2*((x(1,1)*log((x(1,1)/phiie(1,1))))+x(2,1)*log(x(2,1)/phiie(2,1)))+2*x(2,1)*x(1,1)*(((phiie(1,1)/x(1,1))*(1-(lamda(1,1)/lamda(2,1))))+((phiie(2,1)/x(2,1))*(1-(lamda(2,1)/lamda(1,1)))))+x(2,1)*q(1,1)*(((1-((exp(-a(2,1)/temp))*teta(2,1)/((exp(-a(2,1)/temp))*teta(2,1)+teta(1,1)))^2)*log((exp(-a(2,1)/temp)))+(1-(teta(2,1)/((exp(-a(1,1)/temp))*teta(1,1)+teta(2,1)))^2)*(exp(-a(1,1)/temp))*log((exp(-a(1,1)/temp)))))+x(1,1)*q(2,1)*(((1-((exp(-a(1,1)/temp))*teta(1,1)/((exp(-a(1,1)/temp))*teta(1,1)+teta(2,1)))^2)*log(exp(-a(1,1)/temp))+(1-(teta(1,1)/(((exp(-a(2,1)/temp))*teta(2,1))+teta(1,1)))^2)*(exp(-a(2,1)/temp))*log(exp(-a(2,1)/temp)))))-coeffexp);
%x0=[0;0]
x0=[A(2,1);A(3,1)];% affectation au vecteur x0 les valeur initialiser par user 
val=fminsearch( uniff,x0)
a12=val(1,1);
a21=val(2,1);
a=[val(1,1);val(2,1)]
verf1=exp((x(2,1)*log(D0(1,1)))+(x(1,1)*(D0(2,1)))+2*((x(1,1)*log((x(1,1)/phiie(1,1))))+x(2,1)*log(x(2,1)/phiie(2,1)))+2*x(2,1)*x(1,1)*(((phiie(1,1)/x(1,1))*(1-(lamda(1,1)/lamda(2,1))))+((phiie(2,1)/x(2,1))*(1-(lamda(2,1)/lamda(1,1)))))+x(2,1)*q(1,1)*(((1-((exp(-a(2,1)/temp))*teta(2,1)/((exp(-a(2,1)/temp))*teta(2,1)+teta(1,1)))^2)*log((exp(-a(2,1)/temp)))+(1-(teta(2,1)/((exp(-a(1,1)/temp))*teta(1,1)+teta(2,1)))^2)*(exp(-a(1,1)/temp))*log((exp(-a(1,1)/temp)))))+x(1,1)*q(2,1)*(((1-((exp(-a(1,1)/temp))*teta(1,1)/((exp(-a(1,1)/temp))*teta(1,1)+teta(2,1)))^2)*log(exp(-a(1,1)/temp))+(1-(teta(1,1)/(((exp(-a(2,1)/temp))*teta(2,1))+teta(1,1)))^2)*(exp(-a(2,1)/temp))*log(exp(-a(2,1)/temp)))));
erreurr0=(abs(verf1-coeffexp)/coeffexp)*100;vect0=[a12,a21,verf1];
disp('valeur touver par le code du fminsearch')
disp('     a12                     a21                        coeff_trouver                     ')
disp(vect0)
disp('l erreur')
disp(erreurr0)
%% 	section de calcul du 1 er coeff en utilisant les valeurs initialisées par l’utilisateur
%definiton des variable qu on est besoin pour calculer les coefficeints de diffusion possible de different varaition du a12 a21 
k=1;ind=0;inds=0;coefdiffth=coef00;
A(1,1)=sa;%affectation au Vecteur A les valeur initialise par l utilisateur sa = a12 sb= a21
A(2,1)=sb;
A(3,1)=sa;
A(4,1)=sb;
A(5,1)=sa;
A(6,1)=sb;
A(7,1)=sa;
A(8,1)=sb;
delta=1;j=0;

%% 	code spécifique cas ou du coeff exper inferieur du celle calculer en première itération ===>pour avoir convergence

   if coefdiffth > coeffexp
        while coefdiffth-coeffexp ~= 0 & delta == 1
            %code specifique  du variation de a12 et a21 en 4 variations ;err represente le pas 
         
            A(1,1)=A(1,1)-err;%le cas a12 -err
            A(2,1)=A(2,1)-err;%le cas a21 -err
            A(3,1)=A(3,1)+err;%le cas a12 +err
            A(4,1)=A(4,1)+err;%le cas a21 +err
            A(5,1)=A(5,1)+err;%le cas a12 +err
            A(6,1)=A(6,1)-err;%le cas a21 -err
            A(7,1)=A(7,1)-err;%le cas a12 -err
            A(8,1)=A(8,1)+err;%le cas a21 +err
            k=k+1;
            %code pour calcule le coeff theorique  pour chaque cas de variation de a12 a21
            %
                    for i=1:8
                            if mod(i,2)~= 0 
                                %calcule des parametres en fonction de a12 a21 varier
                                %
                            abto(i,1)=exp(-A(i,1)/temp);
                            bato(i,1)=exp(-A(i+1,1)/temp);
                            abteta(i,1)=teta(1,1)*abto(i,1)/(teta(1,1)*abto(i,1)+teta(2,1));
                            bateta(i,1)=teta(2,1)*bato(i,1)/(teta(2,1)*bato(i,1)+teta(1,1));
                            aateta(i,1)=teta(1,1)/(teta(1,1)+teta(2,1)*bato(i,1));
                            bbteta(i,1)=teta(2,1)/(teta(2,1)+teta(1,1)*abto(i,1));
                            part(1,1)=x(2,1)*q(1,1)*(((1-bateta(i,1)^2)*log(bato(i,1))+(1-bbteta(i,1)^2)*abto(i,1)*log(abto(i,1))));
                            part(2,1)=x(1,1)*q(2,1)*(((1-abteta(i,1)^2)*log(abto(i,1))+(1-aateta(i,1)^2)*bato(i,1)*log(bato(i,1))));
                            cas(i,1)=exp(part(1,1)+part(2,1)+part1+part2);
                            inc(k,i)=cas(i,1);
                            indab(k,i)=A(i,1);
                            indba(k,i)=A(i+1,1);
                            
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
                    %detection  valeur plus proche a coef exp par determination du minimum des coeffcient theo calculer par les differnet variations du a12 a21 possible 
                     minimum(k,1)=min(coed);
%                             coefdiffth=minimum(k,1);
                                 %verification que le minimum est superieure coefexp 
                                 %
                                
                                     if minimum(k,1) > coeffexp
                                         coefdiffth=minimum(k,1);
                                         j=j+1;
                                         delta = 1;
                                     else
                                         delta=0;
                                     end
                                 
            

        end
        %% 	identification du ligne et colonne exacte du a12 a21 dans le tableau de variation possible
         for i = 1:k
      for e=1:7
      if coefdiffth==inc(i,e)
          lin=i+1;
          col=e;
      end
      end
  end
   end
   %% 	code spécifique cas ou le coefficient expérimentale est supérieure de celle calculer en première itération le même algorithme et raisonnement pour le cas ou le coefficient expérimentale est supérieure de coef theo
 if coefdiffth < coeffexp
        while coefdiffth-coeffexp ~= 0 & delta == 1
            % variation possible  de a12 a21 
            A(1,1)=A(1,1)-err;
            A(2,1)=A(2,1)-err;
            A(3,1)=A(3,1)+err;
            A(4,1)=A(4,1)+err;
            A(5,1)=A(5,1)+err;
            A(6,1)=A(6,1)-err;
            A(7,1)=A(7,1)-err;
            A(8,1)=A(8,1)+err;
            k=k+1;
            %calcule du parametre to ; teta en fonction des a12 a21 ainsi le coeff de diffusion 
            %
            % 
                    for i=1:8
                            if mod(i,2)~= 0 
                            abto(i,1)=exp(-A(i,1)/temp);
                            bato(i,1)=exp(-A(i+1,1)/temp);
                            abteta(i,1)=teta(1,1)*abto(i,1)/(teta(1,1)*abto(i,1)+teta(2,1));
                            bateta(i,1)=teta(2,1)*bato(i,1)/(teta(2,1)*bato(i,1)+teta(1,1));
                            aateta(i,1)=teta(1,1)/(teta(1,1)+teta(2,1)*bato(i,1));
                            bbteta(i,1)=teta(2,1)/(teta(2,1)+teta(1,1)*abto(i,1));
                            part(1,1)=x(2,1)*q(1,1)*(((1-bateta(i,1)^2)*log(bato(i,1))+(1-bbteta(i,1)^2)*abto(i,1)*log(abto(i,1))));
                            part(2,1)=x(1,1)*q(2,1)*(((1-abteta(i,1)^2)*log(abto(i,1))+(1-aateta(i,1)^2)*bato(i,1)*log(bato(i,1))));
                            cas(i,1)=exp(part(1,1)+part(2,1)+part1+part2);
                            inc(k,i)=cas(i,1);%tableau du du coef calculer pour chaque variation possible de a12 a21
                            indab(k,i)=A(i,1);%tableau du variation a12 
                            indba(k,i)=A(i+1,1);%tableau du variation a21 
                            
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

                           %verification de coeff calcul dapres le variation possible du a12 a21 est proche du coef exp      
                                     if maximum(k,1) < coeffexp
                                         coefdiffth=maximum(k,1);
                                         jj=jj+1;
                                         delta = 1;
                                     else
                                         delta=0;
                                     end
                                 
            

        end
        %% 	identification du ligne et colonne exacte du a12 a21 dans le tableau de variation possible
        for i = 1:k
      for e=1:7
      if coefdiffth==inc(i,e)
          lin=i;
          col=e;
      end
      end
  end
    end
%% localisation et identification du valeur du a12 a21 
  Aab=indab(lin,col);
  Aba=indba(lin,col);
  

   %% vérification par calcule du coeff en utilisant le a12 a21 trouver 
        abto(8,1)=exp(-Aab/temp);
        bato(8,1)=exp(-Aba/temp);
        abteta(8,1)=teta(1,1)*abto(8,1)/(teta(1,1)*abto(8,1)+teta(2,1));
        bateta(8,1)=teta(2,1)*bato(8,1)/(teta(2,1)*bato(8,1)+teta(1,1));
        aateta(8,1)=teta(1,1)/(teta(1,1)+teta(2,1)*bato(8,1));
        bbteta(8,1)=teta(2,1)/(teta(2,1)+teta(1,1)*abto(8,1));
        part(1,1)=x(2,1)*q(1,1)*(((1-bateta(8,1)^2)*log(bato(8,1))+(1-bbteta(8,1)^2)*abto(8,1)*log(abto(8,1))));
        part(2,1)=x(1,1)*q(2,1)*(((1-abteta(8,1)^2)*log(abto(8,1))+(1-aateta(8,1)^2)*bato(8,1)*log(bato(8,1))));
        verf=exp(part(1,1)+part(2,1)+part1+part2);
        coefdiffthf=verf;
        
        A(2,1)=Aab;A(3,1)=Aba;z=1;
        
        
        %% cette section est conçue pour tracer la courbe en fonction de xa
        for w=0:0.01:0.7
            x(1,1)=w;
            x(2,1)=(1-w);dx=x(1,1);
  




lamda=[r(1,1)^(1/3);r(2,1)^(1/3)];
to=[exp(-A(1,1)/temp);exp(-A(2,1)/temp);exp(-A(3,1)/temp);exp(-A(4,1)/temp)];
%a
teta(1,1)=(x(1,1)*(q(1,1)))/(((x(1,1)*(q(1,1)))+(x(2,1)*(q(2,1)))));
phiie(1,1)=(x(1,1)*(lamda(1,1)))/((x(1,1)*(lamda(1,1)))+(x(2,1)*(lamda(2,1))));

%b
teta(2,1)=(x(2,1)*(q(2,1)))/(((x(1,1)*(q(1,1)))+(x(2,1)*(q(2,1)))));
phiie(2,1)=(x(2,1)*(lamda(2,1)))/((x(1,1)*(lamda(1,1)))+(x(2,1)*(lamda(2,1))));
tetabb=teta(2,1)/(to(2,1)*teta(1,1)+teta(2,1));
%ab
tetaab=to(2,1)*teta(1,1)/(to(2,1)*teta(1,1)+teta(2,1));
tetaaa=teta(1,1)/((to(3,1)*teta(2,1))+teta(1,1));
%ba
tetaba=to(3,1)*teta(2,1)/(to(3,1)*teta(2,1)+teta(1,1));
tetabb=teta(2,1)/(to(2,1)*teta(1,1)+teta(2,1));

part1=(x(2,1)*log(D0(1,1)))+(x(1,1)*(D0(2,1)))+2*((x(1,1)*log((x(1,1)/phiie(1,1))))+x(2,1)*log(x(2,1)/phiie(2,1)));
part2=2*x(2,1)*x(1,1)*(((phiie(1,1)/x(1,1))*(1-(lamda(1,1)/lamda(2,1))))+((phiie(2,1)/x(2,1))*(1-(lamda(2,1)/lamda(1,1)))));
part3=x(2,1)*q(1,1)*(((1-tetaba^2)*log(to(3,1))+(1-tetabb^2)*to(2,1)*log(to(2,1))));
part4=x(1,1)*q(2,1)*(((1-tetaab^2)*log(to(2,1))+(1-tetaaa^2)*to(3,1)*log(to(3,1))));
ess(z,1)=exp(part1+part2+part3+part4);w(z,1)=w;
     
       z=z+1 ;
        
        end
    w=linspace(0,0.7,z-1)   ;
    hold on 
    
     
%% section de calcule d'erreur 
erreurr=(abs(verf-coeffexp)/coeffexp)*100;
format longEng
vect=[Aab,Aba,verf];
disp('valeur touver par le code du boucle while et les variation possible de a12 a21 ')
disp('     a12                     a21                        coeff_trouver                    ')
disp(vect)
disp('l erreur')
disp(erreurr)
        A(2,1)=a(1,1);  ;A(3,1)=a(2,1);z=1;

        ess1=0;
        %% cette section est conçue pour tracer la courbe en fonction de xa
        for w=0:0.01:0.7
            x(1,1)=w;
            x(2,1)=(1-w);dx=x(1,1);
  




lamda=[r(1,1)^(1/3);r(2,1)^(1/3)];
to=[exp(-A(1,1)/temp);exp(-A(2,1)/temp);exp(-A(3,1)/temp);exp(-A(4,1)/temp)];
%a
teta(1,1)=(x(1,1)*(q(1,1)))/(((x(1,1)*(q(1,1)))+(x(2,1)*(q(2,1)))));
phiie(1,1)=(x(1,1)*(lamda(1,1)))/((x(1,1)*(lamda(1,1)))+(x(2,1)*(lamda(2,1))));

%b
teta(2,1)=(x(2,1)*(q(2,1)))/(((x(1,1)*(q(1,1)))+(x(2,1)*(q(2,1)))));
phiie(2,1)=(x(2,1)*(lamda(2,1)))/((x(1,1)*(lamda(1,1)))+(x(2,1)*(lamda(2,1))));
tetabb=teta(2,1)/(to(2,1)*teta(1,1)+teta(2,1));
%ab
tetaab=to(2,1)*teta(1,1)/(to(2,1)*teta(1,1)+teta(2,1));
tetaaa=teta(1,1)/((to(3,1)*teta(2,1))+teta(1,1));
%ba
tetaba=to(3,1)*teta(2,1)/(to(3,1)*teta(2,1)+teta(1,1));
tetabb=teta(2,1)/(to(2,1)*teta(1,1)+teta(2,1));

part1=(x(2,1)*log(D0(1,1)))+(x(1,1)*(D0(2,1)))+2*((x(1,1)*log((x(1,1)/phiie(1,1))))+x(2,1)*log(x(2,1)/phiie(2,1)));
part2=2*x(2,1)*x(1,1)*(((phiie(1,1)/x(1,1))*(1-(lamda(1,1)/lamda(2,1))))+((phiie(2,1)/x(2,1))*(1-(lamda(2,1)/lamda(1,1)))));
part3=x(2,1)*q(1,1)*(((1-tetaba^2)*log(to(3,1))+(1-tetabb^2)*to(2,1)*log(to(2,1))));
part4=x(1,1)*q(2,1)*(((1-tetaab^2)*log(to(2,1))+(1-tetaaa^2)*to(3,1)*log(to(3,1))));
ess1(z,1)=exp(part1+part2+part3+part4);w(z,1)=w;
     
       z=z+1 ;
        
        end
   w=linspace(0,0.7,z-1);
figure
 
 plot(w,ess1) 
 title('courpe en utilisant a12 a21 trouver par fminsearch');
 xlabel('fraction xa');
ylabel('coeff diffusion theo');
figure
    plot(w,ess)  
    title('         "courpe en utilisant a12 a21 trouver par boucle while et variation possible    ')
     xlabel('fraction xa');
ylabel('coeff diffusion theo');
 %% remarque 
 % 1
 %quand j initialise fminsearch par a12 =a21 =1 
 % la fonction trouve une valeur de a12 = 234.096488820211e+000 ; a21=1.25762636149561e+003 
 %le coeff_calculer = 13.3000004938735e-006 ; erreur = 3.71333470241173e-006
 % 2
    %valeur touver par le code du boucle while et les variation possible de a12=a21 =1
    %
    %      a12                     a21                        coeff_trouver                     erreur
    %     928.051753968135e+000    928.051753968135e+000    13.2999999999999e-006    420.332139238976e-015
% 3 
        %quand on initialise fminsearch par a12 =a21 = 928
        % valeur touver par le code du fminsearch
        %      a12                     a21                        coeff_trouver                     erreur
        %     928.159309081465e+000    927.997044338426e+000    13.3000000000000e-006    63.6866877634812e-015
 %comme conclusion c est que la valeur initialiser par l utilisateur influence sur la valeur trouver par le code du fminsearch 
 %par contre le code des boucle while et variation par pas converge vers la valeur de 12 a21 ou l'erreur est tres faible
 % 
 %si on prends par example le 1 et 2 de section remarque on constate que fminsearch a obtenue un erreur d ordre 10^(-6)
 %par contre le code de while et des variations a obtenue un erreur d'ordre 10^(-15) 

%% Conclusion
%finalement on peut conclure que la valeur de coefficient de diffusion
%augmente on augmentant la fraction du solute 
% l'allure du courbe est exponentiele croissante 
%pour chaque valeur de a12 a21 initialiser par l'utilisateur  on aura une
%valeur a12 a21 trouver par les deux algos 
%% N.B 
%ces resultats sont trouver pour une valeur de a12 a21 proche des valuer exacte qui verifie la condition que  coeff thero  proche du coeff exp 
%de point de vue le temps d execution qu il soit reduite  car le pas de
%variation est d ordre 1*10^-12 et qu l 'erreur soit d'ordre 10^-15