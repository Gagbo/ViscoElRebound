function [GTaux,VRest, LTemp] = Mashup(varargin)
%Mashup This function takes as arguments :
%       The list of files to analyse
%       The maximum number of peaks to look for picmax
%       Number of different temperatures nT

%%%%Variables initialisation%%%%
nT = varargin{nargin};
picmax = varargin{nargin-1};
Gtaux = zeros(nargin-2,2.*nT,2);
LTemp = zeros(nT);
LTempraw = zeros(nT);
VRest = zeros(nargin-2,2.*nT,picmax-2,2);
Tauxmoy = zeros(nT);
ErreurRT = zeros(nT);
Sortie = zeros(nT,3);

%%%%Temperature list recuperation%%%%
i=1;
k=1;
while LTempraw(nT)==0
	file = sscanf(varargin{i},'%d %*[-] %d %*[-] %d %[-] %d');
	for k = 1:nT
		if file(2)==LTempraw(k)
			break
		elseif LTempraw(k)==0
			LTempraw(k) = file(2);
			disp(['Temperature n ',num2str(k),' is ',num2str(file(2))]);
			break
		end
	end
	i=i+1;
	if i==nargin-1 && LTempraw(nT)==0
		error('Probleme dans la recup de temperature. A priori, c''est que nT n''est pas assez grand dans l''appel de Mashup');
	end
end

LTemp=sort(LTempraw);

%%%%Concat Matrix Filling%%%%
PremierIndice = ones(nT);
for i=1:nargin-2
	j=0; %j is the index of the temperature of the measure in LTemp

	file = sscanf(varargin{i},'%d %*[-] %d %*[-] %d %[-] %d');
	for k=1:nT
		if file(2)==LTemp(k)
			j=k;
			break
		end
	end
	if j==0
		error('Temperature not found in LTemp. nT may be too small');
	end
	[pic,taux,alpha,vitesse] = Anapicsv2(varargin{i},picmax);
	if taux(2)>0.1
		disp([filename,' not taken into account.']);
		continue;
	end
	GTaux(PremierIndice(j),2*j-1,:) = [taux(1),file(3)];
	GTaux(PremierIndice(j),2*j,:) = [taux(2),file(3)];
	VRest(PremierIndice(j),2*j-1,1:numel(alpha),1)=alpha;
	VRest(PremierIndice(j),2*j,1:numel(vitesse),1)=vitesse;
	VRest(PremierIndice(j),2*j,1:numel(alpha),2)=file(3);
	VRest(PremierIndice(j),2*j,1:numel(vitesse),2)=file(3);

	PremierIndice(j) = PremierIndice(j) + 1;
end

%%%%Post-traitement%%%%

%Courbe Taux_de_restit = f(T)
for j=1:nT
	Tauxmoy(j) = sum(GTaux(:,2*j-1,1),1)/(PremierIndice(j)-1);
	ErreurRT(j) = std(GTaux(1:PremierIndice(j)-1,2*j-1,1))/2;
	Sortie(j,2) = Tauxmoy(j);
	Sortie(j,3) = ErreurRT(j);
end
h = figure;
LTempcorig=zeros(nT);
for i=1:nT
	LTempcorig(i) = LTemp(i)+3.1*10^(-6)*LTemp(i)^3-0.0014*LTemp(i)^2+0.2*LTemp(i)-5.3;
	Sortie(i,1) = LTempcorig(i);
end
C1TA = errorbar(LTempcorig,Tauxmoy,ErreurRT);
hold all;
%C150 = plot([LTempcorig(1),LTempcorig(nT)],[0.5,0.5]);
NbEssais=zeros(nT);
%Possibility to plot all the measures, with colors
%according to each sample used.
for j=1:nT
	NbEssais(j)=PremierIndice(j)-1;
	for i=1:PremierIndice(j)-1
		hold all;
		switch GTaux(i,2*j-1,2)
		case 1
			couleur='r';
		case 2
			couleur='g';
		case 4
			couleur='b';
		case 5
			couleur='y';
		case 6
			couleur='k';
		end
		%scatter(LTempcorig(j),GTaux(i,2*j-1,1),[couleur,'+']);
	end
end
xlabel('Temperature (Celsius)','FontSize',15);
ylabel('Coef. of restitution','FontSize',15);
%title('Coef. of restitution with respect to temperature');
ylim([0 1]);
xlim([-0.5 200]);

print(h,'-depsc','-r300','/run/media/apaloo/Dropbox/Stage M1/Images/Resultats/RestTemp9010.eps');
dlmwrite('/run/media/apaloo/Dropbox/Stage M1/Images/Resultats/Temp-Rest-Std.txt', Sortie, 'delimiter', '\t', 'precision', 4, 'newline', 'unix');

hold off;

disp(['Total measures compiled : ',num2str(sum(NbEssais))]);
for i = 1:nT
	disp(['Measures done at ',num2str(LTempcorig(i)),' degrees : ',num2str(NbEssais(i))]);
end
end
