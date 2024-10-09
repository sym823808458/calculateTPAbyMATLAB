%20211210sym更新用于1119数据集的td任务的读取跃迁偶极矩矩阵任务，任务是第一直接提出跃迁偶极矩矩阵并把特征矩阵写入文件，二是算出双光子吸收
clc ; clear all ; close all ;
delete(gcp('nocreate'));
% parpool("local")
set(0,'defaultfigurecolor','w');
FolderName = uigetdir;
cd(FolderName);
filelist = cat(1,dir('*.txt'));
SaveName3 = ['tdTDMFeature','.txt1'];
fid3 = fopen(SaveName3, 'w');
writeline=['samplename',' ','Somax ',' ','Stmax',' ','Sesmax',' ','MUes',' ','MUo',' ',...
    'MUs1',' ','MUoi',' ','MUif',' ','MUmax',' ','Es1',' ','Eomax',' ','Etmax'];
fprintf(fid3, '%s\n', writeline);
fclose(fid3);
tic
for ii =1:length(filelist)
    close all ;
    fid = fopen(filelist(ii).name,'r');
    file=filelist(ii).name
    samplename=filelist(ii).name(1:end-4);
    line_1 = fgetl(fid);  % read line excluding newline character
    groundstatedipmom=str2num(line_1(40:end-4));
    while ~feof(fid)
        tline = fgetl(fid);
        disp(tline)
    end
    statesnum=str2num(tline);
    statesnum=statesnum(1)+1;
    musquare20=zeros(statesnum,statesnum,length(ii));
    statesnumlim=statesnum;%statesnum;
    
    qe=1.602*10^(-19);%C, charge of electron
    me=9.1094*10^(-31);%kg mass of electron
    Debye2Cm=3.33564*10^(-30);
    epsilon0=8.8542*10^(-12);
    c=2.998*10^8;
    au2Cm=qe*0.528*10^(-10);
    h_bar=1.0546*10^(-34);%Js
    eV2radpers=qe/h_bar;
    % tau=0.1*eV2radpers;%rad/s
    tau=0.1*eV2radpers
    tline=[];
    
    allexcitedenergies=zeros(statesnum,1);
    mu=zeros(statesnum,statesnum,3);
    Osscilatorstrength=zeros(statesnum,statesnum);
    i=1;
    mu(1,1,:)=groundstatedipmom;
    stateindex=zeros(statesnum,1);
    stateindex(1)=0;
    fid=fopen(file,'r');
    
    while ~strcmp(tline,' Transition dipole moment between ground state (0) and excited states (a.u.)')
        tline=fgetl(fid);
    end
    tline=fgetl(fid);
    while 1
        tline=fgetl(fid);
        if tline==' '
            fprintf('finished');
            break
        else
            tline=strtrim(tline);
            tnum=str2num(tline);
            i=i+1;
            mu(1,i,:)=tnum(3:5);
            mu(i,1,:)=tnum(3:5);
            Osscilatorstrength(1,i)=tnum(7);
            Osscilatorstrength(i,1)=tnum(7);
            thisenergy=tnum(6);
            stateindex(i)=tnum(1);
            allexcitedenergies(i)=thisenergy;
        end
    end  
    energy=allexcitedenergies*eV2radpers;
    
    while ~strcmp(tline,' Note: In below output the case of i=j corresponds to electronic contribution to dipole moment of excited state i')
        tline=fgetl(fid);
    end
    tline=fgetl(fid);
    tline=fgetl(fid);
    stateindex=0:statesnum-1;
    while 1
        tline=fgetl(fid);
        if tline==-1
            fprintf('finished');
            break
        else
            tline=strtrim(tline);
            tnum=str2num(tline);
            mu(find(stateindex==tnum(1)),find(stateindex==tnum(2)),:)=tnum(3:5);
            mu(find(stateindex==tnum(2)),find(stateindex==tnum(1)),:)=tnum(3:5);
            Osscilatorstrength(find(stateindex==tnum(1)),find(stateindex==tnum(2)))=tnum(7);
            Osscilatorstrength(find(stateindex==tnum(2)),find(stateindex==tnum(1)))=tnum(7);
        end
    end
    fclose(fid);
    
    musquare=mu(:,:,1).^2+mu(:,:,2).^2+mu(:,:,3).^2;

    mu=mu*au2Cm;
mudebye=mu/Debye2Cm;

mudebyesquare=mudebye(:,:,1).^2+mudebye(:,:,2).^2+mudebye(:,:,3).^2;
% musquare=mu(:,:,1).^2+mu(:,:,2).^2+mu(:,:,3).^2;
musquare20=mudebyesquare(1:21,1:21);

%     mudebye=(mu/Debye2Cm).^2;
%     si=size(mudebye)
%     musquare20(:,:,ii)=musquare(1:si(1),1:si(1));
%     musquare20(:,statesnum+1,ii)=allexcitedenergies;
%     muexcitedstate=zeros(statesnum,1);
%     MUo=musquare20(1,1,ii);
%     for j =1:statesnum  %对角元扣掉，之前需要算个读出MUo
%         muexcitedstate=musquare20(j,j,ii);
%         musquare20(j,j,ii)=0;
%     end
%     
%     %写一个有关于mu^2的偶极矩矩阵，单位是au2
%     [MUes,Sesmax] = max(muexcitedstate);
%     [MUoi,Somax] = max(musquare20(1,2:end-1,ii));
%     [MUif,Stmax] = max(musquare20(Somax+1,Somax+2:end-1,ii));
%     Stmax=Stmax+Somax;
%     Eomax=musquare20(Somax+1,end,ii);
%     Etmax=musquare20(Stmax+1,end,ii);
%     MUs1=musquare20(1,2,ii);
%     Es1=musquare20(2,end,ii);
%     [MUmax,S]=max(max(musquare20(2:end,2:end-1,ii)));
%     
%     fid3 = fopen(SaveName3, 'a+');
%     writeline=['samplename',' ','Somax ',' ','Stmax',' ','Sesmax',' ','MUes',' ','MUo',' ',...
%         'MUs1',' ','MUoi',' ','MUif',' ','MUmax',' ','Es1',' ','Eomax',' ','Etmax'];
%     writeline=[samplename,' ',num2str(Somax),' ',num2str(Stmax),' ',num2str(Sesmax),' ',num2str(MUes),' ',num2str(MUo),...
%         ' ',num2str(MUs1),' ',num2str(MUoi),' ',num2str(MUif),' ',num2str(MUmax),' ',num2str(Es1),' ',num2str(Eomax),' ',num2str(Etmax)];
%     fprintf(fid3, '%s\n', writeline);
%     fclose(fid3);
%     
    f=figure;
    yvalues = {'S0','S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20'};
    xvalues = yvalues;
    h = heatmap(xvalues,yvalues,musquare20,'Colormap',parula,'FontName','Arial','FontSize',12)
    titlee= ['heatmap of transition dipole moment'];
    h.Title=titlee;
    h.YLabel='States number';
    h.XLabel='States number';
    formatSpec = ['TDM of ',samplename];
    saveas(f,formatSpec,'jpeg');
  
    %% plot linear absorption spectrum
    % tau=0.5*10^15;%rad/s

    omegaplot=(energy(2)-2*tau)/10^15:0.001:(energy(statesnumlim)+2*tau)/10^15;
    omegaplot=[omegaplot energy(2:statesnumlim)'/10^15];
    omegaplot=sort(omegaplot);
    
    sigma1st=zeros(length(omegaplot),1);
    alpha=zeros(length(omegaplot),3,3);
    
    for i=1:length(omegaplot)
        for x1=1:3
            for x2=1:3
                for j=2:statesnumlim
                    alpha(i,x1,x2)= alpha(i,x1,x2) + 1/h_bar*pi*mu(1,j,x1).*mu(1,j,x2)*(1/(energy(j)-omegaplot(i)*10^15-1i*tau)+1/(energy(j)+omegaplot(i)*10^15+1i*tau));
                end
            end
        end
    end
    plotossicf=[energy(2:end)/10^15 Osscilatorstrength(2:end,1)];
    sigma1st=4*pi^2/c/epsilon0*omegaplot'*10^15.*squeeze(imag(alpha(:,1,1)+alpha(:,2,2)+alpha(:,3,3))/3);
    %convert to wavelength
    [omega1st A1 A2]=unique(omegaplot);
    sigma1stl=sigma1st(A1);
    lambda1=linspace(c/10^6/(omegaplot(1)/2/pi),c/10^6/(omegaplot(end)/2/pi),1000);
    epsilon=interp1(c/10^6./(omega1st/2/pi),sigma1stl*10^4/(3.825*10^(-19)),lambda1);
    %% calculating gama
    for i=1:size(mu,1)
        mu(i,i,:)=mu(i,i,:)-mu(1,1,:);
    end
    % omega=2*pi*299.8/760*10^15; %radpersecond
    omega=[((energy(2)-2*tau)/10^15/2:0.5*tau/10^15:(energy(statesnumlim)+2*tau)/10^15)/2 energy(2:statesnumlim)'/2/10^15];
    % omega=[((energy(2)-2*tau)/10^15/2:0.05:(energy(statesnumlim)+2*tau)/10^15)/2 energy(2:statesnumlim)'/2/10^15];
    %sym changes 2*tau/10^15 to 0.005 for more points
    %omegaplot=(energy(2)-2*tau)/10^15:0.001:(energy(statesnumlim)+2*tau)/10^15;
    omega=sort(omega);
    omega(omega>(energy(2)-2*tau)/10^15)=[];
    omega1set=omega*10^15;
    omega2set=omega*10^15;
    omega3set=-omega*10^15;
    omegasigma=omega1set+omega2set+omega3set;
    
    gama1=zeros(length(omega),3,3,3,3);
    %用于MATLAB并行计算，任务小不需要的时候注释掉这一行把下一行的parfor改成for
    for ry=1:length(omega)
        for x1=1:3
            for x2set=1:3
                for x3set=1:3
                    for x4set=1:3
                        for permutaion=1:6
                            if permutaion==1; x2=x2set; x3=x3set; x4=x4set; omega1=omega1set; omega2=omega2set; omega3=omega3set; end
                            if permutaion==2; x2=x3set; x3=x2set; x4=x4set; omega1=omega2set; omega2=omega1set; omega3=omega3set; end
                            if permutaion==3; x2=x2set; x3=x4set; x4=x3set; omega1=omega1set; omega2=omega3set; omega3=omega2set; end
                            if permutaion==4; x2=x4set; x3=x3set; x4=x2set; omega1=omega3set; omega2=omega2set; omega3=omega1set; end
                            if permutaion==5; x2=x3set; x3=x4set; x4=x2set; omega1=omega2set; omega2=omega3set; omega3=omega1set; end
                            if permutaion==6; x2=x4set; x3=x2set; x4=x3set; omega1=omega3set; omega2=omega1set; omega3=omega2set; end
                            for i=2:statesnumlim
                                for j=2:statesnumlim
                                    for k=2:statesnumlim
                                        gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) + mu(1,i,x1)*mu(i,j,x4)*mu(j,k,x3)*mu(k,1,x2)/((energy(i)-omegasigma(ry)-1i*tau)*(energy(j)-omega1(ry)-omega2(ry)-1i*tau)*(energy(k)-omega1(ry)-1i*tau));
                                        gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) + mu(1,i,x4)*mu(i,j,x1)*mu(j,k,x3)*mu(k,1,x2)/((energy(i)+omega3(ry)+1i*tau)*(energy(j)-omega1(ry)-omega2(ry)-1i*tau)*(energy(k)-omega1(ry)-1i*tau));
                                        %%20211020wave发现code中“(energy(i)-omega3(ry)+1i*tau)”写错，应该是(energy(i)+omega3(ry)+1i*tau)
                                        gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) + mu(1,i,x2)*mu(i,j,x3)*mu(j,k,x1)*mu(k,1,x4)/((energy(i)+omega1(ry)+1i*tau)*(energy(j)+omega1(ry)+omega2(ry)+1i*tau)*(energy(k)-omega3(ry)-1i*tau));
                                        gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) + mu(1,i,x2)*mu(i,j,x3)*mu(j,k,x4)*mu(k,1,x1)/((energy(i)+omega1(ry)+1i*tau)*(energy(j)+omega1(ry)+omega2(ry)+1i*tau)*(energy(k)+omegasigma(ry)+1i*tau));
                                    end
                                    gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) - mu(1,i,x1)*mu(i,1,x4)*mu(1,j,x3)*mu(j,1,x2)/((energy(i)-omegasigma(ry)-1i*tau)*(energy(i)-omega3(ry)-1i*tau)*(energy(j)-omega1(ry)-1i*tau));
                                    gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) - mu(1,i,x1)*mu(i,1,x4)*mu(1,j,x3)*mu(j,1,x2)/((energy(i)-omega3(ry)-1i*tau)*(energy(j)+omega2(ry)+1i*tau)*(energy(j)-omega1(ry)-1i*tau));
                                    gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) - mu(1,i,x4)*mu(i,1,x1)*mu(1,j,x2)*mu(j,1,x3)/((energy(i)+omegasigma(ry)+1i*tau)*(energy(i)+omega3(ry)+1i*tau)*(energy(j)+omega1(ry)+1i*tau));
                                    gama1(ry,x1,x2set,x3set,x4set)= gama1(ry,x1,x2set,x3set,x4set) - mu(1,i,x4)*mu(i,1,x1)*mu(1,j,x2)*mu(j,1,x3)/((energy(i)+omega3(ry)+1i*tau)*(energy(j)-omega2(ry)-1i*tau)*(energy(j)+omega1(ry)+1i*tau));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    gama1=gama1/6/h_bar^3;
    gama=zeros(length(omega),1);
    for i=1:3
        for j=1:3
            if i~=j
                gama = gama + 1/15*squeeze(gama1(:,i,i,j,j)+gama1(:,i,j,i,j)+gama1(:,i,j,j,i));
            end
        end
    end
    for i=1:3
        gama = gama + 1/5*squeeze(gama1(:,i,i,i,i));
    end
    sigma=3/2*h_bar/epsilon0^2*(omega'*10^15).^2/c^2.*imag(gama)/10^(-58);
    lambda2=c/10^6./(omega/2/pi);
    sigmal=sigma;
    xwzy1=2*pi*c./(1e15.*plotossicf(:,1))*1e9';
    ywzy1=plotossicf(:,2);
       
    %plot together with linear absorption
    F=figure
    yyaxis left
    set(gca,'FontName','Times New Roman','linewidth',1.2,'fontsize',12)
    plot(lambda1,epsilon,'linewidth',1.8);
    hold on
    stem(xwzy1,ywzy1*(max(epsilon)/max(ywzy1)),'k-','linewidth',1.8);
    ylabel('molar extinction coefficient/M^{-1}cm^{-1}','FontWeight','bold','FontSize',15)
    yyaxis right
    plot(lambda2/2,sigmal,'r','linewidth',1.8);
    ylabel('TPA cross section/sigma GM','FontWeight','bold','FontSize',15);
    xlabel('wavelength/nm','FontWeight','bold','FontSize',15','FontWeight','bold');
    
    % xlim([min(lambda1) max(lambda2/2)]);
    ylim([0 max(sigmal)*1.4]);
    xlim([min(lambda1) max(lambda1)]);
    title(strcat('TPA and linear absorption spectra of ',samplename));
    plotfororiginx1=lambda2/2;
    plotfororiginxy2=ywzy1*(max(epsilon)/max(ywzy1));
    
    dd=allexcitedenergies;
    dd(dd==0) = [];
    dd=1240./dd;
    B=[plotfororiginx1' sigmal];
    TPAplot=zeros(length(dd),3);
    TPAplot(:,1)=linspace(1,length(dd),length(dd)) ;
    TPAplot(:,2)=dd;
    for p =1:length(dd)
        f1=dd(p)-plotfororiginx1;
        [max2,Index2] = min(abs(f1));
        TPAplot(p,3)=sigmal(Index2);
    end
    hold on
    plot(TPAplot(:,2),TPAplot(:,3),'r*','linewidth',1.8);
    legend('linear absorption','Oscillator Strength','TPA plot with \lambda/2','states','location','northwest');
    legend('boxoff')
    fx = diff(TPAplot(:,3));  
    for pp=1:length(fx)
        if TPAplot(pp,3)==TPAplot(pp+1,3) %jianbingtai
            disp(['S',num2str(pp),'and S',num2str(pp+1),'is degenerated states'])
        end
        if pp<=length(fx)-3
            if fx(pp)>0 & fx(pp+1)>0 &fx(pp+2)<0&fx(pp+3)<0
                str=['S',num2str(pp+2)];
                text(TPAplot(pp+2,2),TPAplot(pp+2,3),str);
                disp(['The two photon state is ',str])
            end
        end
    end

    [dx,dy]=max(sigmal);
    [~,Index] = min(abs(dd-plotfororiginx1(dy)));%Index返回索引
    formatSpec = [samplename,' S',num2str(Index)];
    str=['S',num2str(Index)];
    text(dd(Index),dx,str);
    saveas(F,formatSpec,'jpeg');
    fprintf(['The max two photon state is ',num2str(Index)])
    filename = ['TPAcalforplot',samplename,'.xlsx'];
    A = [lambda1' epsilon'];
    Nr = normalize(epsilon','range');
    osi_cor=max(Osscilatorstrength(1,:))*Nr;
    A = [lambda1' osi_cor];
    sheet = 1;
    xlswrite(filename,A,sheet);
    xlRange = 'C1';
    xlswrite(filename,B,sheet,xlRange);
    C=[xwzy1 plotfororiginxy2/max(plotfororiginxy2)];
    xlRange = 'E1';
    xlswrite(filename,C,sheet,xlRange);
    hold off
    clear F
end
toc