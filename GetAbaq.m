function [Output] = GetAbaq(ID,varargin)
%GetAbaq This function takes as arguments :
%	Identification of output files
%	The list of files to analyse

%%%%Variables initialisation%%%%
nargin
nT = floor(nargin/2)
Output = zeros(nT,5);
LTemp = zeros(nT,1);
LTempRaw = zeros(nT,1);

%%%%Temperature list recuperation%%%%
i=1;
k=1;
while LTempRaw(nT)==0
	Cond = sscanf(varargin{i},'%c %*[-] %d %*[-] %d');
	for k = 1:nT
		if Cond(2)==LTempRaw(k)
			break
		elseif LTempRaw(k)==0
			LTempRaw(k) = Cond(2);
			disp(['Temperature n°',num2str(k),' is ',num2str(Cond(2))]);
			break
		end
	end
	i=i+1;
	if i == nargin-1 && LTempRaw(nT)==0
		error('Temperature recuperation');
	end
end

LTemp = sort(LTempRaw)

%%%%Concat Matrix Filling%%%%
for i = 1:nargin-1
	Cond = sscanf(varargin{i},'%c %*[-] %d %*[-] %d');
	j = 0;
	for k =1:nT
		if Cond(2)==LTemp(k)
			j=k;
			break
		end
	end
	if j==0
		error('Temperature not found in LTemp');
	end
	if Cond(1)=='E'
		[Dura,Peak,Int] = AnaAbaqEff(varargin{i});
		Output(j,2) = Peak;
		Output(j,3) = Dura;
		Output(j,4) = Int;
	elseif Cond(1)=='S'
		[COR] = AnaAbaqSpeed(varargin{i});
		Output(j,5) = COR;
	end
end

%%%%Post treatment%%%%
for i = 1:nT
	Output(i,1) = LTemp(i);
end

h=figure;
subplot(2,2,1);
Curve1 = plot(Output(:,1),Output(:,2));
hold all;
title([ID,' Force peak of impact with respect to temperature']);
xlabel('Temperature (Celsius)');
ylabel('Peak (N)');
hold off;

subplot(2,2,2);
Curve2 = plot(Output(:,1),Output(:,3));
hold all;
title([ID,' Duration of impact with respect to temperature']);
xlabel('Temperature (Celsius)');
ylabel('Duration of impact (s)');
hold off;

subplot(2,2,3);
Curve3 = plot(Output(:,1),Output(:,4));
hold all;
title([ID, ' Area under the impact with respect to temperature']);
xlabel('Temperature (Celsius)');
ylabel('Area (N.s)');
hold off;

subplot(2,2,4);
Curve4 = plot(Output(:,1),Output(:,5));
hold all;
title([ID,' COR with respect to temperature']);
xlabel('Temperature (Celsius)');
ylabel('COR');
CORexp = dlmread('/home/gerry/Dropbox/Stage M1/Images/Resultats/Temp-Rest-Std.txt');
errorbar(CORexp(:,1),CORexp(:,2),CORexp(:,3));
hold off;

print(h,'-dpng','-r400',['/home/gerry/Dropbox/Stage M1/Images/Resultats/Efforts/Courbes-',ID,'.png']);

dlmwrite(['/home/gerry/Dropbox/Stage M1/Images/Resultats/Efforts/Temp-Peak-Dura-Area-COR-',ID,'.txt'],[Output(:,1),Output(:,2),Output(:,3),Output(:,4),Output(:,5)], 'Delimiter', '\t', 'Precision', 4, 'newline', 'unix');

end

function [Dura,Peak,Int] = AnaAbaqEff(filename)

%%%%Variable initialisation%%%%
Peak = 0;
Dura = 0;
BegI = 1;
EndI = 1;
Int = 0;

%%%%File reading%%%%
FileID = fopen(filename);
File = textscan(FileID,'%f %f %*[^\n]','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',3,'BufSize',10000000);
fclose(FileID);
File = [File{1},File{2}];
Siz = length(File);

%%%%Data analysis%%%%
Peak = File(1,2);
PeakInd = 1;
StartFound = 0;
EndFound = 0;

for i = 1:Siz
	if File(i,2) > Peak
		PeakInd = i;
		Peak = File(i,2);
	end

	if ~StartFound && File(i,2) > 0
		BegI = i;
		StartFound = 1;
		disp('Start of impact found');
	end

	if StartFound && ~EndFound && File(i,2) == 0
		EndI = i;
		EndFound = 1;
		disp('End of impact found');
	end

end

if ~StartFound
	disp('Start of impact not found');
elseif ~EndFound
	disp('End of impact not found');
end

Dura = File(EndI,1) - File(BegI,1);
Int = trapz(File(BegI:EndI,1),File(BegI:EndI,2));

disp([filename,' : the peak is at ',num2str(Peak),' N']);
disp([filename,' : the duration of impact is ',num2str(Dura),' s']);
disp([filename,' : the area under impact is ',num2str(Int),' N.s']);

end


function [COR] = AnaAbaqSpeed(filename)

%%%%File reading%%%%
FileID = fopen(filename);
File = textscan(FileID,'%f %f %*[^\n]','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',3,'BufSize',10000000);
fclose(FileID);
File = [File{1},File{2}];
Siz = length(File);

COR = abs(File(Siz,2)/File(1,2));

end

