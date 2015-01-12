function [Output] = Mashimpact(varargin)
%Mashimpact This function takes as arguments :
%	The list of files to analyse
%	Identification of output files
%	Number of different temperatures nT

%%%%Variables initialisation%%%%
nT = varargin{nargin};
ID = varargin{nargin-1};
Data = zeros(nT,3,10);
LTemp = zeros(nT);
LTempRaw = zeros(nT);
Output = zeros(nT,4,2);

%%%%Temperature list recuperation%%%%
i=1;
k=1;
while LTempRaw(nT)==0
	Cond = sscanf(varargin{i},'%d %*[-] %d %*[-] %d %*[-] %d');
	for k = 1:nT
		if Cond(2)==LTempRaw(k)
			break
		elseif LTempRaw(k)==0
			LTempRaw(k) = Cond(2);
			disp(['Temperature n ',num2str(k),' is ',num2str(Cond(2))]);
			break
		end
	end
	i=i+1;
	if i == nargin-1 && LTempraw(nT)==0
		error('Temperature recuperation : surely nT not big enough');
	end
end

LTemp = sort(LTempRaw);

%%%%Concat Matrix Filling%%%%
FirstIndex = ones(nT);
for i = 1:nargin-2
	[NewFile,Dura,Peak,Int,Cond] = Anaimpact(varargin{i});
	j=0;
	for k=1:nT
		if Cond(2)==LTemp(k)
			j=k;
			break
		end
	end
	if j==0
		error('Temperature not found in LTemp. nT may be too small');
	end
	Data(j,1,FirstIndex(j)) = Peak;
	Data(j,2,FirstIndex(j)) = Dura;
	Data(j,3,FirstIndex(j)) = Int;
	FirstIndex(j) = FirstIndex(j) + 1;
end

%%%%Post treatment%%%%

for i = 1:nT
	for j=1:3
		Output(i,j+1,1) = sum(Data(i,j,:),3)/(FirstIndex(i)-1);
		Output(i,j+1,2) = std(Data(i,j,1:FirstIndex(i)))/2;
	end
	Output(i,1,1) = LTemp(i)+3.1*10^(-6)*LTemp(i)^3-0.0014*LTemp(i)^2+0.2*LTemp(i)-5.3;
end

Outsim = dlmread('/path/to/Temp-Peak-Dura-Area-COR-Sim69-4n.txt');
Sizsim = size(Outsim,1);

h1 = figure;
%subplot(2,2,1);
Sim1 = plot(Outsim(:,1),Outsim(:,2)./Outsim(Sizsim,2));
hold all;
Curve1 = errorbar(Output(:,1,1),Output(:,2,1)./Output(nT,2,1),Output(:,2,2));
%title(['Normalized peak of impact with respect to temperature']);
xlabel('Temperature (Celsius)','FontSize',15);
ylabel('Peak','FontSize',15);

hold off;

h2 = figure;
%subplot(2,2,2);
Sim2 = plot(Outsim(:,1),Outsim(:,3));
hold all;
Curve2 = errorbar(Output(:,1,1),Output(:,3,1),Output(:,3,2));
%title(['Duration of impact with respect to temperature']);
xlabel('Temperature (Celsius)','FontSize',15);
ylabel('Duration of impact (s)','FontSize',15);
hold off;

h3 = figure;
%subplot(2,2,3);
Sim3 = plot(Outsim(:,1),Outsim(:,4)./Outsim(Sizsim,4));
hold all;
Curve3 = errorbar(Output(:,1,1),Output(:,4,1)./Output(nT,4,1),Output(:,4,2));
%title(['Area under the impact, Hot value-normalized, with respect to temperature']);
xlabel('Temperature','FontSize',15);
ylabel('Area','FontSize',15);
hold off;

h4 = figure;
%subplot(2,2,4);
Curve4 = plot(Outsim(:,1),Outsim(:,5));
hold all;
%title([ID,' COR with respect to temperature']);
xlabel('Temperature (Celsius)','FontSize',15);
ylabel('COR','FontSize',15);
CORexp = dlmread('/path/to/Temp-Rest-Std.txt');
errorbar(CORexp(:,1),CORexp(:,2),CORexp(:,3));
xlim([0 200]);
hold off;


print(h1,'-depsc','-r300',['/path/to/',ID,'-Peak.eps']);
print(h2,'-depsc','-r300',['/path/to/',ID,'-Dura.eps']);
print(h3,'-depsc','-r300',['/path/to/',ID,'-Area.eps']);
print(h4,'-depsc','-r300',['/path/to/',ID,'-COR.eps']);

%dlmwrite(['/home/gerry/Dropbox/Stage M1/Images/Resultats/Efforts/Temp-Peak-Std-',ID,'.txt'], [Output(:,1,1),Output(:,2,1),Output(:,2,2)], 'Delimiter', '\t', 'Precision', 4, 'newline', 'unix');
%dlmwrite(['/home/gerry/Dropbox/Stage M1/Images/Resultats/Efforts/Temp-Dura-Std-',ID,'.txt'], [Output(:,1,1),Output(:,3,1),Output(:,3,2)], 'Delimiter', '\t', 'Precision', 4, 'newline', 'unix');
%dlmwrite(['/home/gerry/Dropbox/Stage M1/Images/Resultats/Efforts/Temp-Int-Std-',ID,'.txt'], [Output(:,1,1),Output(:,4,1),Output(:,4,2)], 'Delimiter', '\t', 'Precision', 4, 'newline', 'unix');
end
