function z = Sheep_wolves_grass(nfile)
%
if nargin ~=1
    nfile='';
else
    nfile=num2str(nfile); % Simulation number to be used in the output files
end

%% Parameters
par.sheep0 =     100; % Initial number of Sheep
par.wolves0 =     50; % initial number of Wolves
par.grassgrowth = 30; % Number of timesteps required for a patch to turn green
par.sheepgain =    4; % Energy gain by Sheep for eating one patch of grass
par.sheepgrowth =  4; % Probability of Sheep generating offspring
par.wolfgain =    20; % Energy gain by Wolves for eating one Sheep
par.wolfgrowth =   5; % Probability of Wolves generating offspring
par.worldx =      51; % x and y dimensions of the world
par.worldy =      51;
par.maxsheep =  30000; % Used to pre-allocate tables, but not required. Although it will run slower if exceeded during the simulation
par.maxwolves = 30000;
par.maxtime =     501; % Number of timesteps  0 .. maxtime-1
par.graf =          1; % Graphic output? 1- Yes  2- No
par.addnum =        1; % Number of simulations averaged

load('GSWcolormap'); % color scheme for the output

for k = 1:par.addnum  % k number of simulations being averaged
    
    
    %% Setup grass
    %Grass is initiated either green (0) or brown (1), and randomly given a
    %number between 1 and par.grassgrowth to set the clock of the brown
    %grass. At each timestep the clock is reduced by 1 until a patch of grass
    %becomes available to eaten (green).
    grass = randi([0 1], par.worldx, par.worldy) .* floor(1 + rand(par.worldx, par.worldy) * (par.grassgrowth - 1));
        
    %% create sheep table
    %                              1        2    3 4    5      6
    %                            alive, heading, x,y,energy, patch
    
    sheep = nan(par.maxsheep,6); % Sheep table pre allocation
    sheep(1:par.sheep0,1) = 1;   % The first par.sheep0 are set to alive
    sheep(1:par.sheep0,2) = rand(par.sheep0,1) * 2*pi();                 % given an initial random heading
    sheep(1:par.sheep0,3) = rand(par.sheep0,1) * par.worldx;             % given a random x location
    sheep(1:par.sheep0,4) = rand(par.sheep0,1) * par.worldy;             % given a random y location
    sheep(1:par.sheep0,5) = rand(par.sheep0,1) * 2 * par.sheepgain;      % given a random energy level
    sheep(:,6) = floor(sheep(:,3)) + 1 + floor(sheep(:,4)) * par.worldx; % based on x,y assigned to a patch number
    
    %% create wolves table
    %                               1        2    3 4    5      6
    %                             alive, heading, x,y,energy, patch
    
    wolves = nan(par.maxwolves,6); % Wolves table pre allocation
    wolves(1:par.wolves0,1) = 1;   % The first par.wolves0 are set to alive
    wolves(1:par.wolves0,2) = rand(par.wolves0,1) * 2*pi();                 % given an initial random heading
    wolves(1:par.wolves0,3) = rand(par.wolves0,1) * par.worldx;             % given a random x location
    wolves(1:par.wolves0,4) = rand(par.wolves0,1) * par.worldy;             % given a random y location
    wolves(1:par.wolves0,5) = rand(par.wolves0,1) * 2 * par.wolfgain;       % given a random energy level
    wolves(:,6) = floor(wolves(:,3)) + 1 + floor(wolves(:,4)) * par.worldx; % based on x,y assigned to a patch number
    
    time=0; %initial time step
    if isfinite(par.maxtime)     % Pre alocation of the output table if maxtime is not infinite
        output = nan(par.maxtime,4); 
    else
        output=[];
    end
    
    while time<par.maxtime % timesteps for each simulation
        
        % Generating output for current timestep
        green = (grass == 0); % Determines which patches are green
        
        if isfinite(par.maxtime) % Add current state of the simulation to the output

            %    Time, Total green grass, total live sheep, total live wolves
            output(time+1,:) = [time, sum(sum(green)), sum(sheep(:,1) > 0), sum(wolves(:,1) > 0)]; 
        else
            output = [output; time, sum(sum(green)),sum(sheep(:,1) > 0),sum(wolves(:,1) > 0)];
        end
                
        %% Plotting
        if par.graf == 1
            figure(1)
            
            livesheep =   (sheep(:,1)>0);
            livewolves = (wolves(:,1)>0);
            
            mapa = (grass(:)+1);
            mapa(sheep(livesheep,6)) =   35;
            mapa(wolves(livewolves,6)) = 36;
            mapa=reshape(mapa,par.worldx,par.worldy);
            mapa(par.worldx,par.worldy) = 256;
            mapa = transpose(mapa);
            
            heatmap(mapa)
            colormap(GSWCM)
            colorbar('off')
            title('Wolves - red; Sheep - white; Grass - green; Growing grass - dark shades of green')
            
            figure(2) % Sheep energy histogram
            histogram(sheep(sheep(:,1)>0,5),20)
            title(strcat('Histogram of sheep energy. Total energy=',num2str(sum(sheep(:,1)>0))))

            figure(3) % Wolves energy histogram
            histogram(wolves(wolves(:,1)>0,5),20)
            title(strcat('Histogram of wolves energy. Total energy=',num2str(sum(wolves(:,1)>0))))

            figure(4) % Grass clock histogram
            histogram(grass,30)
            title(strcat('Histogram of grass growth. Total green grass=',num2str(sum(sum(green)))))

            figure(5) % Current output
            plot(output(:,1),output(:,2:end))
            xlabel('Simulation time')
            ylabel('Number of individuals')
            legend('Grass','Sheep','Wolves')
        end
        
        %% Sheep
        if sum(sheep(:,1)>0)>0 % if there are live sheep
            
            % move
            sheep(:,2) = sheep(:,2) + rand(size(sheep,1),1) * (50 / 360 * 2*pi()) - rand(size(sheep,1),1) * (50 / 360 * 2*pi());
            sheep(:,3) = sheep(:,3) + cos(sheep(:,2));
            sheep(:,4) = sheep(:,4) + sin(sheep(:,2));
            
            sheep((sheep(:,3) < 0),3) =          sheep((sheep(:,3) < 0         ),3) + par.worldx;
            sheep((sheep(:,3) > par.worldx),3) = sheep((sheep(:,3) > par.worldx),3) - par.worldx;
            
            sheep((sheep(:,4) < 0),4) =          sheep((sheep(:,4) < 0         ),4) + par.worldy;
            sheep((sheep(:,4) > par.worldy),4) = sheep((sheep(:,4) > par.worldy),4) - par.worldy;
            
            sheep(:,6) = floor(sheep(:,3)) + 1 + floor(sheep(:,4))*par.worldx;
            
            %lose energy
            sheep(sheep(:,1) > 0,5) = sheep(sheep(:,1) > 0,5) - 1;
            
            %eat-grass
            tobefeed = find(sheep(:,1) > 0);
            sequence = randperm(length(tobefeed)); %randomize feeding sequence
            grasstemp = grass(:);
            for m = 1 : length(tobefeed) % Goes through each sheep
                mm = sequence(m);
                if  grasstemp(sheep(tobefeed(mm),6)) == 0  %checks if patch is green
                    grasstemp(sheep(tobefeed(mm),6)) = par.grassgrowth; %grass turns brown
                    sheep(tobefeed(mm),5) = sheep(tobefeed(mm),5) + par.sheepgain; %sheep gains energy
                end
            end
            grass = reshape(grasstemp, par.worldx, par.worldy);
            
            
            %sheep death
            deadsheep = find(sheep(:,1) > 0 & sheep(:,5) < 0); %Lists live sheep with negative energy
            if ~isempty(deadsheep)
                for m = 1 : length(deadsheep)
                    sheep(deadsheep(m),1) = 0; %sets to dead
                    sheep(deadsheep(m), 2:6) = nan; %and blanks all properties
                end
            end
            
            % reproduce-sheep
            livesheep = find(sheep(:,1) > 0);
            for m = 1 : length(livesheep)
                if rand * 100 < par.sheepgrowth
                    emptysheep = find(sheep(:,1) ~= 1);
                    sheep(livesheep(m), 5) = sheep(livesheep(m), 5) / 2; %half the parent energy
                    sheep(emptysheep(1),:) = sheep(livesheep(m),:); %The parent gets cloned into the offspring
                    
                    %Offspring gets a random direction and moves
                    sheep(emptysheep(1),2) = rand*2*pi(); 
                    sheep(emptysheep(1),3) = sheep(emptysheep(1),3) + cos(sheep(emptysheep(1),2)); 
                    sheep(emptysheep(1),4) = sheep(emptysheep(1),4) + sin(sheep(emptysheep(1),2));
                end %if
            end% for
            
            sheep((sheep(:,3) < 0         ),3) = sheep((sheep(:,3) < 0         ),3) + par.worldx;
            sheep((sheep(:,3) > par.worldx),3) = sheep((sheep(:,3) > par.worldx),3) - par.worldx;
            
            sheep((sheep(:,4) < 0         ),4) = sheep((sheep(:,4) < 0         ),4) + par.worldy;
            sheep((sheep(:,4) > par.worldy),4) = sheep((sheep(:,4) > par.worldy),4) - par.worldy;
            
            sheep(:,6) = floor(sheep(:,3)) + 1 + floor(sheep(:,4)) * par.worldx;
        end
        %% Wolves
        
        if sum(wolves(:,1) > 0) > 0  % if there are live wolves
            % move
            wolves(:,2) = wolves(:,2) + rand(size(wolves,1),1) * (50 / 360 * 2*pi()) - rand(size(wolves,1),1) * (50 / 360 * 2*pi());
            wolves(:,3) = wolves(:,3) + cos(wolves(:,2));
            wolves(:,4) = wolves(:,4) + sin(wolves(:,2));
            
            wolves((wolves(:,3) < 0         ),3) = wolves((wolves(:,3) < 0         ),3) + par.worldx;
            wolves((wolves(:,3) > par.worldx),3) = wolves((wolves(:,3) > par.worldx),3) - par.worldx;
            
            wolves((wolves(:,4) < 0         ),4) = wolves((wolves(:,4) < 0         ),4) + par.worldy;
            wolves((wolves(:,4) > par.worldy),4) = wolves((wolves(:,4) > par.worldy),4) - par.worldy;
            
            wolves(:,6) = floor(wolves(:,3)) + 1 + floor(wolves(:,4)) * par.worldx;
            
            %lose energy
            wolves(wolves(:,1) > 0,5) = wolves(wolves(:,1) > 0,5) - 1;
            
            %eat sheep
            livewolves = find(wolves(:,1) > 0);  %List live wolves
            sequence = randperm(length(livewolves)); %random sequence for feeding
            for m = 1 : length(livewolves)  %For each wolve that needs to feed
                mm = sequence(m);
                sheephere = find(sheep(:,6) == wolves(mm,6)); %Lists all sheep with the same possition as the wolf
                if ~isempty(sheephere) %if there are sheep near the wolf
                    seque2 = randperm(length(sheephere));
                    wolves(mm,5) = wolves(mm,5) + par.wolfgain; %wolf gain energy
                    sheep(sheephere(seque2(1)),1) = 0; %Sheep dies
                    sheep(sheephere(seque2(1)),2:6) = nan; %Sheep properties are blanked
                end %if
            end %for
            
            %wolves death
            deadwolves = find(wolves(:,1) > 0 & wolves(:,5) < 0);
            if ~isempty(deadwolves)
                for m = 1 : length(deadwolves)
                    wolves(deadwolves(m),1) = 0;
                    wolves(deadwolves(m),2:6) = nan;
                end
            end
            
            % reproduce-wolves
            livewolves = find(wolves(:,1) > 0);
            for m = 1 : length(livewolves)
                if rand * 100 < par.wolfgrowth
                    emptywolves = find(wolves(:,1) ~= 1);
                    wolves(livewolves(m),5) =  wolves(livewolves(m),5) / 2; %half the parent energy
                    wolves(emptywolves(1),:) = wolves(livewolves(m),:); %offspring=parent
                    
                    wolves(emptywolves(1),2) = rand * 2*pi(); %move offspring
                    wolves(emptywolves(1),3) = wolves(emptywolves(1),3) + cos(wolves(emptywolves(1),2));
                    wolves(emptywolves(1),4) = wolves(emptywolves(1),4) + sin(wolves(emptywolves(1),2));
                end %if
            end %for
            
            wolves((wolves(:,3) < 0         ),3) = wolves((wolves(:,3) < 0         ),3) + par.worldx;
            wolves((wolves(:,3) > par.worldx),3) = wolves((wolves(:,3) > par.worldx),3) - par.worldx;
            
            wolves((wolves(:,4) < 0         ),4) = wolves((wolves(:,4) < 0         ),4) + par.worldy;
            wolves((wolves(:,4) > par.worldy),4) = wolves((wolves(:,4) > par.worldy),4) - par.worldy;
            
            wolves(:,6) = floor(wolves(:,3)) + 1 + floor(wolves(:,4)) * par.worldx;
        end
        %% Grass
        
        %grow grass
        grass = grass - 1; %decreases the clock on all patches
        grass(grass < 0) = 0; 
        
        if sum(wolves(:,1) > 0) + sum(sheep(:,1) > 0) < 1
            if sum(green(:)) >= par.worldx*par.worldy
                break % break simulation if there are no animals and grass is all green
            end
        end
        time = time + 1; % add 1 to the timestep
    end %while time<par.maxtime
    
    if par.addnum > 1
        clear deadsheep
        clear deadwolves
        clear emptysheep
        clear emptywolves
        clear livesheep
        clear livewolves
        clear m
        clear mapa
        clear mm
        clear sequence
        clear seque2
        clear sheephere
        clear tobefeed
        clear grass
        clear green
        clear sheep
        clear wolves
        clear grasstemp
    end
    
    if exist('Noutput','var')
        Noutput=Noutput+output; %adds output to the previous simulations
    else
        Noutput=output; 
    end
    
end % for k (next simulation if adding several)

%%
%Generate final output and record files
z = Noutput/k; % calculates the mean of the outputs of the k simulations

Gn = z(:,2);
Sn = z(:,3);
Wn = z(:,4);

OutputFolder = '';
save(strcat(OutputFolder,'Grass',nfile,'.txt'),'Gn','-ascii')
save(strcat(OutputFolder,'Sheep',nfile,'.txt'),'Sn','-ascii')
save(strcat(OutputFolder,'Wolves',nfile,'.txt'),'Wn','-ascii')

end % function