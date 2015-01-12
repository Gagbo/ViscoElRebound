function [File,Mastercurve,CWLF,TMes,Shi] = GetWLFv2(filename,Tref)

%%%%Global variables%%%%
global MEST; %Measures done per temperature in the DMA
MEST = 16;
global NbMe;

%%%%File opening and other variables initialisation%%%%
File = dlmread([filename,'.txt'], '\t', 1,1);

[NbMe,NbVa] = size(File);
Freq = File(:,1); % Frequency of the DMA measure
Temp = File(:,2); % Temperature of the DMA measure
Gpri = File(:,3); % Storage Modulus
Gsec = File(:,4); % Loss Modulus
Delt = File(:,5); % Damping Factor
Torq = File(:,6); % Torque applied

Mastercurve = zeros(NbMe,7);
% Columns
% Frequency associated with Tref for the measure
% Temperature
% Storage Modulus
% Damping Factor
% Frenquency-shift with respect to Tref measures
% base 10 logarithm of the Shift

%%%%Tref indice search%%%%
Iref = 0;
i = 1;
while (i<NbMe) && (Iref==0)
	if Temp(i) == Tref
		Iref = i
		Irefmod = (MEST - 1 + i)/MEST
	end
	i = i+1;
end

%%%%Mastercurve Initialisation%%%%

for i = 0:MEST-1
	for j = 1:5
		Mastercurve(Iref+i,j) = File(Iref+i,j);
	end
	Mastercurve(Iref+i,6) = 1;
	Mastercurve(Iref+i,7) = 0;
end
disp(['Mastercurve initialised']);

%%%%Right shitfting%%%%

Imod = 0;
while Imod < Irefmod-1
	disp(['Right shift n° ',num2str(Imod+1),' on ',num2str(Irefmod-1),' started']);
	StRef = (Irefmod-1-Imod)*MEST + 1;
	[Shft,Datashft]=ShiftCurve(Mastercurve(StRef:StRef+MEST-1,:),File(StRef-MEST:StRef-1,:));
	Mastercurve(StRef-MEST:StRef-1,:) = Datashft;
	Mastercurve(StRef-MEST:StRef-1,6) = Mastercurve(StRef:StRef+MEST-1,6).*Mastercurve(StRef-MEST:StRef-1,6);
	Mastercurve(StRef-MEST:StRef-1,7) = log(Mastercurve(StRef-MEST:StRef-1,6));
	Imod = Imod + 1;
end

disp(['Right shift ended']);

%%%%Left Shitfing%%%%

Imod = 0;
while Imod < NbMe/MEST - Irefmod
	disp(['Left shift n° ',num2str(Imod+1),' on ',num2str(NbMe/MEST - Irefmod),' started']);
	StRef = (Irefmod - 1 + Imod)*MEST + 1;
	[Shft,Datashft]=ShiftCurve(Mastercurve(StRef:StRef+MEST-1,:),File(StRef+MEST:StRef+2*MEST-1,:));
	Mastercurve(StRef+MEST:StRef+2*MEST-1,:) = Datashft;
	Mastercurve(StRef+MEST:StRef+2*MEST-1,6) = Mastercurve(StRef+MEST:StRef+2*MEST-1,6).*Mastercurve(StRef:StRef+MEST-1,6);
	Mastercurve(StRef+MEST:StRef+2*MEST-1,7) = log(Mastercurve(StRef+MEST:StRef+2*MEST-1,6));
	Imod = Imod + 1;
end

disp(['Left shift ended']);

%%%%Post-treatment%%%%

i=1;
j=1;
TMes = zeros([NbMe/MEST-1,1]);
WLFY = zeros([NbMe/MEST-1,1]);
WLFX = zeros([NbMe/MEST-1,1]);
LShi = zeros([NbMe/MEST-1,1]);
while i < NbMe + 1
	if File(i,2) == Tref
		i = i + MEST;
		continue;
	else
		TMes(j) = File(i,2);
		i = i + MEST;
		j = j + 1;
	end
end

j=1;
for i = 1:NbMe/MEST-1
	WLFX(i) = 1./(TMes(i)-Tref);
end
for i = 1:NbMe/MEST
	if i == Irefmod
		continue;
	else
		WLFY(j) = 1./Mastercurve((i-1)*MEST+1,7);
		LShi(j) = Mastercurve((i-1)*MEST+1,7);
		Shi(j) = Mastercurve((i-1)*MEST+1,6);
		j = j+1;
	end
end

h = figure;
subplot(2,2,1);
%Curve = scatter(Mastercurve(:,1),Mastercurve(:,3),'r+');
set(gca, 'XScale', 'log','YScale', 'log');
hold all;

Mastercurve=sortrows(Mastercurve);

Curve2 = plot(Mastercurve(:,1),Mastercurve(:,3));
set(gca, 'XScale', 'log','YScale', 'log');
title(['Storage modulus (Pa) with respect to frequency at Tref = ',num2str(Tref)]);


subplot(2,2,2);
Curve3 = semilogx(Mastercurve(:,1),Mastercurve(:,5));
title(['tan(delta) with respect to frequency at Tref = ',num2str(Tref)]);


subplot(2,2,3);
TMes(NbMe/MEST) = Tref;
LShi(NbMe/MEST) = 0;
Shi(NbMe/MEST) = 1;
WLFX
WLFY
[FitL,S] = polyfit(WLFX,WLFY,1)
CWLF(1) = -1./FitL(2);
CWLF(2) = FitL(1)./FitL(2);
funcWLF = @(T) (CWLF(1)*(Tref-T)./(CWLF(2)+(T-Tref)));
X = Temp(1):.1:Temp(NbMe);
WLF = funcWLF(X);
Curve4 = plot(X,WLF);
hold all;
scatter(TMes,LShi,'r+');
ylim([-50 50]);
xlim([25 95]);
title('Verification WLF Law');
xlabel('Temperature');
ylabel('log(shift)');
hold off;

subplot(2,2,4);
scatter(WLFX,WLFY,'g+');
hold all;
FitY = FitL(1).*WLFX+FitL(2);
Curve5 = plot(WLFX,FitY);
title('Linear fit for WLF verification');
xlabel('1/(T-Tref)');
ylabel('-1/log(shift)');
text('Units','normalized','Position',[0.5 0.25],'String',char('y = A * x + B', ['A = ',num2str(FitL(1))],['B = ',num2str(FitL(2))]));
print(h,'-dpng','-r1000',['/home/gerry/Dropbox/Stage M1/Images/Resultats/',filename,'/',filename,'-',num2str(Tref),'.png']);
hold off;

[Tabdel,FreqEqu] = DeltaShift(Mastercurve,CWLF,88,Tref,TMes,Shi);

end

function [Shft,Datashft]=ShiftCurve(DataRef,DataMob)

global MEST;

Datashft = zeros(MEST,7);
Targ = floor(MEST/2);

%DataRef(:,1)
%DataRef(:,3)

Ref = pchip(DataRef(:,1),DataRef(:,3));
Mob = pchip(DataMob(:,1),DataMob(:,3));

Obj = DichotSolve(DataRef,DataMob(Targ,3),[DataRef(1,1),DataRef(MEST,1)]);
disp(['Dichotomy done']);
Shft = Obj - DataMob(Targ,1); %The shift is just the one to match the mid-range frequency point on the master curve already done.

%%%%Filling of the Mastercurve w/ modified frequencies%%%%
Datashft(Targ,1) = DataMob(Targ,1)+Shft; %New frequency
Datashft(Targ,2) = DataMob(Targ,2); %Temperature, Gpri, Gsec and damping factor untouched
Datashft(Targ,3) = DataMob(Targ,3);
Datashft(Targ,4) = DataMob(Targ,4);
Datashft(Targ,5) = DataMob(Targ,5);
Datashft(Targ,6) = Datashft(Targ,1)/DataMob(Targ,1); %Shift-factor
Datashft(Targ,7) = log(Datashft(Targ,6));

disp(['log(Shift-factor) = ',num2str(Datashft(Targ,7))]);


for i = [1:Targ-1,Targ+1:MEST]
	Datashft(i,1) = DataMob(i,1)*Datashft(Targ,6);
	Datashft(i,2) = DataMob(i,2);
	Datashft(i,3) = DataMob(i,3);
	Datashft(i,4) = DataMob(i,4);
	Datashft(i,5) = DataMob(i,5);
	Datashft(i,6) = Datashft(Targ,6);
	Datashft(i,7) = log(Datashft(i,6));
end

end

function [Answ]=DichotSolve(Data,Const,Domain)

%This function solves f(X) = Const w/ regard to X, where f is represented by Data(:,1)
%and Data(:,2). Domain is the first domain which will be cut by Dichotomy.

Func = pchip(Data(:,1),Data(:,3));
Doma = Domain;

Answ = mean(Doma);
Rech = (Doma(2)-Doma(1))/2;

while abs(ppval(Func,Answ)-Const) > 0.001*Const
	if abs(ppval(Func,Answ+Rech)-Const)>abs(ppval(Func,Answ-Rech)-Const)
		Answ = Answ-Rech/2;
	else
		Answ = Answ+Rech/2;
	end
	Rech = Rech/2;
	if abs(ppval(Func,Answ+Rech)-ppval(Func,Answ-Rech)) < 0.0001 * Const
		Rech = Rech * 4;
	end
	%disp([num2str(ppval(Func,Answ)-Const), ' ' , num2str(Const)]);
	%disp([num2str(Answ),' ',num2str(Rech)]);
end
end

