% Matlab script to plot the results from the Julia model output
warning off

% Make figures directory
mkdir figures

% Get output file
fname = ['../out/out_addW_20240213.nc'];

% Get array variables to load
vars = {'p','b','z','n','d','o','pIC','bIC','zIC','nIC','dIC','oIC'};

% Get grid settings
grid_vars = {'H','dz'};

% Load into struct
for i = 1:length(vars)
    out.(vars{i}) = ncread(fname,vars{i});
end
for i = 1:length(grid_vars)
    out.(grid_vars{i}) = ncread(fname,grid_vars{i});
end

% Extract individual bacteria
bact = {'aer','nar','nai','nao','nir','nio','nos','aoa','nob','aox'};
for i = 1:length(bact)
    out.(bact{i}) = squeeze(out.b(:,i,:));
    out.([bact{i},'IC']) = squeeze(out.bIC(:,i));
end

% Extract individual N tracers
ncycle = {'nh4','no2','no3','n2o','n2'};
for i = 1:length(ncycle)
    out.(ncycle{i}) = squeeze(out.n(:,i,:));
    out.([ncycle{i},'IC']) = squeeze(out.nIC(:,i));
end

% Extract O2
out.o2   = squeeze(out.o);
out.o2IC = squeeze(out.oIC);

% Get grid cell centers
zc = ((out.dz./2):out.dz:(out.H - (out.dz./2)));

% Get nrec
nrec = size(out.p,3);

% Get colormap for old timesteps
old_cmap = cmocean('balance',nrec+30);
old_cmap = old_cmap(16:end-15,:);

% N + O2 SPINUP
if (1)
    % Get vars,titles
    cmd  = 'rm figures/N_spinup*png';
    system(cmd)
    vars = {'o2','nh4','no2','no3','n2o','n2'};
    tits = {'O$_2$','NH$^{+}_4$','NO$^{-}_2$','NO$^{-}_3$','N$_2$O','N$_2$'};
   
    % Initialize figure
    fig    = piofigs('lfig',0.75);
    sb(1)  = subplot(3,7,1);
    sb(2)  = subplot(3,7,2);
    sb(3)  = subplot(3,7,3);
    sb(4)  = subplot(3,7,4);
    sb(5)  = subplot(3,7,5);
    sb(6)  = subplot(3,7,6);

    % Shrink position, get limits
    for i = 1:length(vars)
        sb(i).Position(4) = sb(i).Position(4)*0.95;
        pos{i} = sb(i).Position;
        tmp = out.(vars{i});
        lims{i} = [0 max(tmp(:))];
    end
    
    % Make initial conditions
    for i = 1:length(vars)
        set(fig,'CurrentAxes',sb(i));
        tmp = out.([vars{i},'IC']);
        plot(tmp,zc,'-','linewidth',1,'color',rgb('Black'));
        hold on
        xlim([lims{i}]);
        ylim([0 2000]);
        set(gca,'YDir','Reverse');
        title(tits{i},'Interpreter','Latex');
        if i > 1
            set(gca,'YTickLabel',[]);
        end
        set(gca,'FontSize',6);
        grid on
        sb(i).Position = pos{i};
    end
    export_fig('-png',['figures/N_spinup_00'],'-m5');
    close all

    % Now go through all timesteps
    for t = 1:nrec
        disp([' ']);
        disp(['tstep = ',num2str(t)]);
        disp([' ']);
        % Initialize figure
        fig    = piofigs('lfig',0.75);
        sb(1)  = subplot(3,7,1);
        sb(2)  = subplot(3,7,2);
        sb(3)  = subplot(3,7,3);
        sb(4)  = subplot(3,7,4);
        sb(5)  = subplot(3,7,5);
        sb(6)  = subplot(3,7,6);
        for i = 1:length(vars)
            set(fig,'CurrentAxes',sb(i));
            for tt = 1:t
                plot(out.(vars{i})(:,tt),zc,'-','linewidth',0.5,'color',old_cmap(tt,:));
                hold on
            end
            plot(out.(vars{i})(:,t),zc,'-','linewidth',1,'color',rgb('Black'));
            hold on
            xlim([lims{i}]);
            ylim([0 2000]);
            set(gca,'YDir','Reverse');
            title(tits{i},'Interpreter','Latex');
            if i > 1
                set(gca,'YTickLabel',[]);
            end
            set(gca,'FontSize',6);
            grid on
            sb(i).Position = pos{i};
        end
        if t < 10
            export_fig('-png',['figures/N_spinup_0',num2str(t)],'-m5');
        else
            export_fig('-png',['figures/N_spinup_',num2str(t)],'-m5');
        end
        close all
    end
end

% BIOMASS SPINUP
if (1)
    % Get title colors
    heter_clr = rgb('Indigo');
    chemo_clr = rgb('FireBrick');
    bio_clr   = rgb('DarkGreen');

    % Get p,z,d, + bacterial names, titles, limits
    vars  = {'p',...
             'z',...
             'd',...
             'aer',...
             'nar',...
             'nai',...
             'nao',...
             'nir',...
             'nio',...
             'nos',...
             'aoa',...
             'nob',...
             'aox'};
    tits1 = {'Phyto C',...
             'Zoo C',...
             'Org C',...
             'AER C',...
             'NAR C',...
             'NAI C',...
             'NAO C',...
             'NIR C',...
             'NIO C',...
             'NOS C',...
             'AOA C',...
             'NOB C',...
             'AOX C'};
    tits2 = {' ',...
             ' ',...
             ' ',...
             ' ',...
             'NO$^{-}_3$ $\rightarrow$ NO$^{-}_2$',...
             'NO$^{-}_3$ $\rightarrow$ N$_2$O' ,...
             'NO$^{-}_3$ $\rightarrow$ N$_2$',...
             'NO$^{-}_2$ $\rightarrow$ N$_2$O',...
             'NO$^{-}_2$ $\rightarrow$ N$_2$',...
             'N$_2$O $\rightarrow$ N$_2$',...
             'NH$^{+}_4$ $\rightarrow$ NO$^{-}_2$',...
             'NO$^{-}_2$ $\rightarrow$ NO$^{-}_3$',...
             'O$^{-}_2$ + NH$^{+}_4$ $\rightarrow$ N$_2$'};
    clrs =   [bio_clr;...
              bio_clr;...
              bio_clr;...
              heter_clr;...
              heter_clr;...
              heter_clr;...
              heter_clr;...
              heter_clr;...
              heter_clr;...
              heter_clr;...
              chemo_clr;...
              chemo_clr;...
              chemo_clr];

    % Get max bacterial biomass (for limits)
    for i = 1:length(vars)
        tmp   = squeeze(out.(vars{i}));
        lims{i} = [0 max(tmp(:))];
    end

    % Initialize figure
    fig    = piofigs('lfig',0.75);
    sb(1)  = subplot(3,7,3);
    sb(2)  = subplot(3,7,4);
    sb(3)  = subplot(3,7,5);
    sbcnt  = 3;
    for i = 1:7
        sb(i+3) = subplot(3,7,7+i);
        sbcnt   = sbcnt + 1;
    end
    for i = 1:3
        sb(i+10) = subplot(3,7,16+i); 
        sbcnt   = sbcnt + 1;
    end

    % Shrink
    for i = 1:sbcnt
        sb(i).Position(4) = sb(i).Position(4)*0.95;
        pos{i} = sb(i).Position;
    end

    % Get initial conditions
    for i = 1:length(vars)
        set(fig,'CurrentAxes',sb(i));
        plot(out.([vars{i},'IC']),zc,'-','linewidth',1,'color',rgb('Black'));
        hold on
        xlim([lims{i}]);
        ylim([0 2000]);
        set(gca,'YDir','Reverse');
        title({tits1{i},tits2{i}},'Interpreter','Latex','color',clrs(i,:));
        if ~ismember(i,[1 4 11]);
            set(gca,'YTickLabel',[]);
        end 
        set(gca,'FontSize',6);
        grid on
        sb(i).Position = pos{i};
    end
    export_fig('-png',['figures/C_spinup_00'],'-m5');
    close all

    % Now go through all timesteps
    for t = 1:nrec
        disp([' ']);
        disp(['tstep = ',num2str(t)]);
        disp([' ']);
        % Initialize figure
        fig    = piofigs('lfig',0.75);
        sb(1)  = subplot(3,7,3);
        sb(2)  = subplot(3,7,4);
        sb(3)  = subplot(3,7,5);
        for i = 1:7
            sb(i+3) = subplot(3,7,7+i);
        end
        for i = 1:3
            sb(i+10) = subplot(3,7,16+i); 
        end

        % Get data up until tstep 
        for i = 1:length(vars)
            set(fig,'CurrentAxes',sb(i));
            this_var = squeeze(out.(vars{i}));
            for tt = 1:t
                plot(out.(vars{i})(:,tt),zc,'-','linewidth',0.5,'color',old_cmap(tt,:));
                hold on
            end
            plot(out.(vars{i})(:,t),zc,'-','linewidth',1,'color',rgb('Black'));
            xlim([lims{i}]);
            ylim([0 2000]);
            set(gca,'YDir','Reverse');
            title({tits1{i},tits2{i}},'Interpreter','Latex','color',clrs(i,:));
            if ~ismember(i,[1 4 11]);
                set(gca,'YTickLabel',[]);
            end
            set(gca,'FontSize',6);
            grid on
            sb(i).Position = pos{i};
        end
        if t < 10
            export_fig('-png',['figures/C_spinup_0',num2str(t)],'-m5');
        else
            export_fig('-png',['figures/C_spinup_',num2str(t)],'-m5');
        end
        close all
    end
end
