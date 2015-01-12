function [pic,taux,alpha,vitesse] = Anapicsv2(filename,picmax)
%Anapics This function takes the .wav file
%        and gives an array with the number of the peak associated
%        with a time.

addpath(genpath('.'));

file = sscanf(filename,'%d %*[-] %d %*[-] %d %[-] %d');
[niv,fq] = wavread(filename,'double');
m = length(niv); %m is the number of samples
g = 9.81; 
tpsrbd = 0.1; %Estimated time between each rebound (s)
perc = 100000; %Progress indicator
global picmaxeff
picmaxeff=picmax;
%picmax is the number of peak to get at most

%Initialisation
Ampmax = 0;
midpics = false;
j=1;
compt=0;
pic(picmax) = 0;
alpha(picmax)=0;
vitesse(picmax)=0;
%%%%Song analysis%%%%
%Search of biggest amplitude
for i = 2:m
    if abs(niv(i)) > Ampmax
        Ampmax = abs(niv(i));
    end
end
disp(['Amplitude maximale : ',num2str(Ampmax)]);

%Sensibility for the peaks calculus
seuil = 50*(1-exp(-abs(file(2)-85)/40))+15;

i=2;
while j<picmax+1 && i<m+1
    if ((abs(niv(i))/Ampmax > seuil/100) && (midpics == false))
        pic(j)=i;
        j=j+1;
        midpics = true;
        compt=0;
    end
    
    if (abs(niv(i))/Ampmax < seuil/100-0.05)
    %if (abs(niv(i)/Ampmax < 0.05))
        compt = compt+1;
    end
    if compt > tpsrbd*fq
        midpics = false;
    end
    
    i = i+1;

    if floor(i/perc) == i/perc
        disp([num2str(i),'sample on ',num2str(m),' reached']);
    end
end

%Recursive call if the script looked for too many peaks
if pic(picmax) ==0
	if picmax==2
		error(['No peak found in ',filename]);
		return
	else
		disp(['Can''t find ',num2str(picmax),' peaks, looking for one less']);
		picmaxeff=picmaxeff-1;
		[pic,taux,alpha,vitesse]=Anapicsv2(filename,picmaxeff);
	end
end

if picmax==picmaxeff


	%Verification (commented because of Mashup calling Anapics relentlessy)
	%Verifpic = figure;
	%x = 1:m;
	%p = plot(x,niv(x));
	%set(gca,'XTick',pic);
	%set(gca,'XTickLabel',{1:picmax});
	%set(p,'Color','red','LineWidth',1);
    %set(gca,'FontSize',15);
    %ylabel('Relative sound level');
	%title([filename,' Peaks verification']);
    %print(Verifpic,'-depsc','-r300','/run/media/apaloo/Dropbox/Stage M1/Images/Resultats/Anapics/BigPeaksVerif.eps');

	%Post-treatment
	for i = 1:picmax
    		pic(i) = pic(i)/fq; %Conversion into seconds
	end


	if picmax == 2
		vitesse(1)=4.0;
		vitesse(2)=g*(pic(2)-pic(1))/2;
		alpha(1)=vitesse(2)/vitesse(1);
		alpha(2)=[];
	else
		alpha(picmax-2)=0; %Initialisation of restitution table
		vitesse(picmax-1)=0; %Initialisation of speed table
		vitesse(1) = 4.0; %We chose 81.5cm as the constant falling height
		%so the first impact speed is 4.0 m.s^-1
	
		for i = 1:picmax-1
	    		vitesse(i+1)=g*(pic(i+1)-pic(i))/2
	    		alpha(i)=vitesse(i+1)/vitesse(i)
	    		%vitesse(i) is the impact speedassociated with alpha(i)
	    		%So vitesse(2) is the speed just before the 2nd impact etc.
		end
		alpha(picmax)=[];
		alpha(picmax-1)=[];
		vitesse(picmax)=[];
		vitesse(picmax-1)=[];
		%Attempt to see correlation between speed and restitution when visquous
		if not(picmax==2 || picmax==3)
			figure;
			p2 = plot(vitesse,alpha);
			xlabel('Impact speed (m/s)');
			ylabel('Coef. of restit associated');
			title([filename,' COR with respect to impact speed']);
		end
	end		
	taux=[mean(alpha),std(alpha)];

	disp([filename,', Mean of the COR measured on ',num2str(picmax),' peaks is ',num2str(taux(1))]);
	disp(['Standard deviation is ',num2str(taux(2))]);
else
	return
end
end
