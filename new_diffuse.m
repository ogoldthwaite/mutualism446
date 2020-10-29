% CS446 -- Computational Modeling and Simulation II
% Owen Goldthwaite, Gautam Mitra, Lolo Niemiec
% September 30, 2020
% Initial Simulation for Mutualism

%% Simulation Parameters %%%
% Seed the random number generator for testing
rng_set = rng(123456789);

% Time-related variables
dt = 1;             % timestep, increment by days, need to go by hours
simLength = 365;    % length of simulation: 1 year
numIterations = 1 + simLength/dt;
animation_fps = 15;  % Speed of visualization

% Grid dimensions
row_count = 50; % width
col_count = 50; % length

%% Constants %%
EMPTY = 1;
GROWING_PLANT = 2;
PLANT = 3;
POLLINATED_PLANT = 4;
ANIMAL = 5;
POLLINATED_ANIMAL = 6;
POLLEN = 7;

prob_init_plant = 0.01; % initial probability a cell is plant
prob_init_animal = 0.01; % initial probability a cell is animal
prob_animal_death = 0.05; % probability of animal death at each timestep
prob_plant_death = 0.05; % probability plant death at each timestep
prob_plant_growth = 0.1; % probability plant growth at each timestep
prob_pollination = 0.1; % probability animal will pick up/drop pollen
prob_pollen_production = 0.8; % chance a plant will produce pollen
prob_pollen_spread = 0.15; % chance pollen spreads to an adjacent empty cell
prob_animal_reproduce = .1; % chance 2 animals reproduce if conditions are good
prob_random_move = 0.15; % chance an animal will move randomly with no stimuli

init_animal_count = 100; % Number of animals at initial animal spawn points
init_plant_count = 100; % Number of plants at initial plant spawn points
init_pollen_count = 100;     % Initial amount of pollen a pollinated plant has
animal_pollen_diffuse_rate = 0.5; % Amount of animals that will diffuse to pollen tile
animal_empty_diffuse_rate = 0.01; % Amount of animals that will diffuse to empty tile
pollen_diffuse_rate = 0.001; % Amount of pollen that will diffuse

%% Counters for statistics
plant_counter = zeros(1, numIterations); % Keep track of the number of plants
animal_counter = zeros(1, numIterations); % Keep track # of animals
animal_pop_counter = zeros(1, numIterations); % Keep track # of animals

%% Set up grids
% Initialize grid to be all empty
grids = ones(row_count, col_count, numIterations) * EMPTY;

animal_pop_grids = zeros(row_count, col_count, numIterations) * EMPTY;
pollen_conc_grids = zeros(row_count, col_count, numIterations) * EMPTY;

for row = 1:row_count
    for col = 1:col_count
        % Grid initialization, plants first
        plant_chance = rand;
        animal_chance = rand;
        
        % Get random number, if smaller than p, add a plant
        if plant_chance < prob_init_plant
            grids(row, col, 1) = PLANT;
            pollen_conc_grids(row, col, 1) = init_pollen_count;
            
        % Now adding animals if smaller than p
        elseif animal_chance < prob_init_animal
            grids(row, col, 1) = ANIMAL;
            animal_pop_grids(row, col, 1) = init_animal_count;
        end
    end
end

% First value for counters
animal_counter(1) = sum(sum(grids(:,:,1)==ANIMAL)) +...
                              sum(sum(grids(:,:,1)==POLLINATED_ANIMAL));
plant_counter(1) = sum(sum(grids(:,:,1)==PLANT)) + ...
                             sum(sum(grids(:,:,1)==GROWING_PLANT)) + ...
                             sum(sum(grids(:,:,1)==POLLINATED_PLANT));
animal_pop_counter(1) = sum(sum(animal_pop_grids(:,:,1)));


pollen_conc_counter(1) = sum(sum(pollen_conc_grids(:,:,1)));


disp("Grids Initialized");

%% Main Simulation Loop
for frame = 2:numIterations
    
    %% Absorbing boundary condition
    % Create grids that are the size of the grid + 4 on each side
    extended_grid_size = size(grids( : , : , frame-1))+4;
    % Initialize all of the cells as empty
    extended_grid = ones(extended_grid_size) * EMPTY; 
    extended_animal_grid = zeros(extended_grid_size) * EMPTY; 
    extended_pollen_grid = zeros(extended_grid_size) * EMPTY; 

    
    % Set the inside portion of the grids equal to the corresponding values 
    % from the previous timestep (a.k.a the previous frame)
    extended_grid(3:end-2, 3:end-2) = grids(:,:,frame-1);
    extended_animal_grid(3:end-2, 3:end-2) = animal_pop_grids(:,:,frame-1);
    extended_pollen_grid(3:end-2, 3:end-2) = pollen_conc_grids(:,:,frame-1);

    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to original(non-extended) grid
    for row = 3:row_count + 1
        for col = 3:col_count + 1

            % Store current cell
            current_cell = extended_grid(row, col);
            current_animal_pop_cell = extended_animal_grid(row,col);
            current_pollen_conc = extended_pollen_grid(row,col);

            updated_cell = current_cell;
            updated_animal_pop_cell = current_pollen_conc;
            updated_pollen_cell = current_pollen_conc;
            
            % Getting Moore Neighborhood
            north = extended_grid(row - 1, col);
            east  = extended_grid(row, col - 1);
            south = extended_grid(row + 1, col);
            west  = extended_grid(row, col + 1);          
            northeast = extended_grid(row - 1, col - 1);
            southeast  = extended_grid(row + 1, col - 1);
            northwest = extended_grid(row - 1, col + 1);
            southwest  = extended_grid(row + 1, col + 1);
            
            % Hello newbie! Getting Von Neumann neighborhood
            northnorth = extended_grid(row - 2, col);
            easteast =   extended_grid(row, col - 2);
            southsouth = extended_grid(row + 2, col);
            westwest =   extended_grid(row, col + 2);
            northeast2 = extended_grid(row - 2, col - 2);
            southeast2 = extended_grid(row + 2, col - 2);
            northwest2 = extended_grid(row - 2, col + 2);
            southwest2 = extended_grid(row + 2, col + 2);
            northleft =  extended_grid(row - 2, col + 1);
            northright = extended_grid(row - 2, col - 1);
            easttop =    extended_grid(row - 1, col - 2);
            eastbot =    extended_grid(row + 1, col - 2);
            southleft =  extended_grid(row + 2, col + 1);
            southright = extended_grid(row + 2, col - 1); 
            westtop =    extended_grid(row - 1, col + 2);
            westbot =    extended_grid(row + 1, col + 2);        

            % Put all of the neighbors into a list
            neighbors = [north, east, south, west, ...
             northeast, southeast, northwest, southwest];
         
            % Put all of the von neighbors into a list
            von_neighbors = [northnorth, easteast, southsouth, westwest,...
                northeast2, southeast2, northwest2, southwest2, ...
                northleft, northright, easttop, eastbot, southleft, ...
                southright, westtop, westbot];

            % List of Neighbor coordinates
            neighbor_coords = [ row-1 col; row col-1; row+1 col; ...
                row col+1; row-1 col-1; row+1 col-1; row-1 col+1; row+1 col+1 ];   
            
            % List of von neighbor coordinates
            von_neighbor_coords = [row-2 col; row col-2; row+2 col; ...
                                   row col+2; row-2 col-2; row+2 col-2;...
                                   row-2 col+2; row+2 col+2; row-2 col+1;...
                                   row-2 col-1; row-1 col-2; row+1 col-2;...
                                   row+2 col+1; row+2 col-1; row-1 col+2;...
                                   row+1 col+2];

            %% Update cell

            % Getting different counts of neighbors
            empty_count = sum(neighbors == EMPTY);
            pollen_count = sum(neighbors == POLLEN);
            norm_plant_count = sum(neighbors == PLANT);
            poll_plant_count = sum(neighbors == POLLINATED_PLANT);
            grow_plant_count = sum(neighbors == GROWING_PLANT);
            norm_animal_count = sum(neighbors == ANIMAL);
            poll_animal_count = sum(neighbors == POLLINATED_ANIMAL);
    
            % Total animal count
            animal_count = norm_animal_count + poll_animal_count;
            % Total plant count
            % Note: pollen and growing plants not included 
            plant_count = poll_plant_count + norm_plant_count;        
            % Total amount of neighbors
            num_neighbors = animal_count + plant_count + pollen_count;

            % Getting total population of all animals in neighboring cells
            if(animal_count > 0)
                neighbor_animal_pop = 0;
                for val = 1:length(neighbors==ANIMAL)
                    neighbor_animal_pop = neighbor_animal_pop + extended_animal_grid(neighbor_coords(val,1), neighbor_coords(val,2));                 
                    neighbor_animal_pop = round(neighbor_animal_pop);
                end
            end
            
            % Getting total amount of pollen in neighboring plant cells
            if (plant_count > 0)
                plant_n_pollen_conc = 0;
                for val = 1:length(neighbors==PLANT)
                    plant_n_pollen_conc = plant_n_pollen_conc + extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                    plant_n_pollen_conc = round(plant_n_pollen_conc);
                end
            end
            
            % Getting total amt of pollen in neighboring pollinated animals
            if (poll_animal_count > 0)
                neighbor_poll_animal_conc = 0;
                for val = 1:length(neighbors==POLLINATED_ANIMAL)
                    neighbor_poll_animal_conc = neighbor_poll_animal_conc + + extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                    neighbor_poll_animal_conc = round(neighbor_poll_animal_conc);
                end
            end        
            
            % Getting total number of all pollen in neighboring cells
            neighbor_pollen_conc = 0;
            for val = 1:length(neighbors)
                neighbor_pollen_conc = neighbor_pollen_conc + extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                neighbor_pollen_conc = round(neighbor_pollen_conc);
            end
            
            % Pollen Diffusion Behavior
            if(current_pollen_conc > 0 || plant_count > 0 || neighbor_pollen_conc > 0)
                % Plants produce pollen and add it to this cells current pollen
                concentration_gain1 = 0;
                if(plant_count > 0)
                    % All neighboring plants have a chance to produce pollen
                    for val = 1:plant_count
                        if(rand < prob_pollen_production)
                            concentration_gain1 = concentration_gain1 + init_pollen_count;
                        end
                    end
                end

                % If there are empty cells than concentration will be decreased by
                % the empty cell count times the diffuse rate since each empty cell
                % will take current concentration * diffuse rate out
                % Only happens if there is no neighboring pollen
                concentration_loss1 = 0;
                if(empty_count > 0 && neighbor_pollen_conc <= 0)
                    concentration_loss1 = current_pollen_conc * ((pollen_diffuse_rate) * empty_count);
                    % If losing more concentration than is possible, just lose current amount
                    % Similar for negative values
                    if(concentration_loss1 > current_pollen_conc)
                        concentration_loss1 = current_pollen_conc;
                    elseif(concentration_loss1 < 0)
                        concentration_loss1 = 0;
                    end
                end

                % If there are pollen neighbors then some pollen will diffuse to
                % lower population neighbors. So either gain or lose concentration
                % based on neighbor populations
                concentration_loss2 = 0;
                concentration_gain2 = 0;
                if(neighbor_pollen_conc > 0)
                    for val = 1:length(neighbors)
                        % Only do this stuff to pollen neighbors
                        if(extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2)) > 0) 
                            cur_neighbor_conc = extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                            % Adding to this cell if neighbor has more pollen
                            if(cur_neighbor_conc > current_pollen_conc)
                                % Gaining concentration, so gain is increased
                                concentration_gain2 = concentration_gain2 + (cur_neighbor_conc * pollen_diffuse_rate);
                            % If current population is greater than neighbor than lose some animals
                            elseif(cur_neighbor_conc < current_pollen_conc)
                                % Losing concentration so concentration loss is increased
                                concentration_loss2 = concentration_loss2 + (current_pollen_conc * pollen_diffuse_rate);
                            end
                        end
                    end
                    if(concentration_loss2 > current_pollen_conc)
                        concentration_loss2 = current_pollen_conc;
                    elseif(concentration_loss2 < 0)
                        concentration_loss2 = 0;
                    end
                end
                % Setting the new populations based on population loss / gain
                updated_pollen_cell = (current_pollen_conc - concentration_loss1 - concentration_loss2 + concentration_gain1 + concentration_gain2);
            end
            
            % Empty cell hello!
            % Don't forget to put checks these are greater than 1!
            if (current_cell == EMPTY)
                if (plant_count > 0 && rand < prob_plant_growth)
                    % Growing plants will commence
                    updated_cell = GROWING_PLANT;
                elseif (plant_count > 0 && rand < prob_pollen_production)
                    % Le plants would like to disperse their pollen now
                    % We should dispere uniformly by using rates, empty now
                    % gets all neighboring plants' pollen*diffuse_rate
                    updated_cell = POLLEN;
                    %updated_pollen_cell = plant_n_pollen_conc * pollen_diffuse_rate;
                elseif (pollen_count > 0 && rand < prob_pollen_spread)
                    % Le pollen will now spread to this empty cell 
                    updated_cell = POLLEN;
                    %updated_pollen_cell = pollen_n_pollen_conc * pollen_diffuse_rate;
                elseif (animal_count > 0 && rand < animal_empty_diffuse_rate)
                    % Animals are boring
                    % They diffuse
                    updated_cell = ANIMAL;
                    updated_animal_pop_cell = neighbor_animal_pop * animal_empty_diffuse_rate;
                elseif (poll_plant_count > 0 && rand < prob_pollen_production)
                    % The pollinated animals would like to get rid of at
                    % least SOME of their pollen can they live?
                    updated_cell = POLLEN;
                    %updated_pollen_cell = poll_plant_n_conc * pollen_diffuse_rate;
                else
                    updated_cell = EMPTY;
                    udpated_animal_pop_cell = 0;
                    %updated_pollen_cell = 0;
                end
                
            % Now for the animal cell
            % Remember to put checks in that enough animals inside!
            elseif (current_cell == ANIMAL)
                if (empty_count > 0 && rand < animal_empty_diffuse_rate)
                    % The animal will diffuse to the empty cell at a rate,
                    % so it loses an amt of pollen, this number based on
                    % how many empty cells are near
                    updated_animal_pop_cell = current_animal_pop_cell - empty_count * animal_empty_diffuse_rate * current_animal_pop_cell;
                    updated_cell = ANIMAL;
                elseif (pollen_count > 0)
                    % There's pollen lurking near, so animal wants to move
                    % there, nothing to do except empty the cell
                    updated_cell = EMPTY;
                    updated_animal_pop_cell = 0;
                elseif (plant_count > 0 && rand < prob_pollination)
                    % If animal is near a plant, then it picks up that good
                    % good pollen, but how much (??) is the begging question
                    updated_cell = POLLINATED_ANIMAL;
                    updated_animal_pop_cell = current_animal_pop_cell;
                    %updated_pollen_cell = init_pollen_count;
                elseif (rand < prob_animal_death)
                    % Rip our dear friend they must go
                    updated_cell = EMPTY;
                    updated_animal_pop_cell = 0;
                else
                    updated_cell = ANIMAL;
                    updated_animal_pop_cell = current_animal_pop_cell;
                end
                
            % Pollinated animal behavior is up and I could care less
            elseif (current_cell == POLLINATED_ANIMAL)
                if (plant_count > 0 && rand < prob_pollination)
                    % Soo if there's a plant nearby, we are going to drop
                    % our pollen off and go back to a normal animal
                    updated_cell = ANIMAL;
                    %updated_pollen_cell = 0;
                    updated_animal_pop_cell = current_animal_pop_cell;
                elseif (pollen_count > 0)
                    % If there is pollen nearby, we will chase it like tom
                    % chases jerry. or jerry chases tom. I can't remember
                    updated_cell = EMPTY;
                    updated_animal_pop_cell = 0;
                    %updated_pollen_cell = 0;
                else
                    updated_cell = POLLINATED_ANIMAL;
                    updated_animal_pop_cell = current_animal_pop_cell;
                    %updated_pollen_cell = current_pollen_conc;
                end
                
            % This is the good stuff
            elseif (current_cell == POLLEN)
                % % We know when the pollen is leaving, either an empty cell
                % % is lurking or we want to diffuse to other pollen cells.
                % % Pollen is increasing when it picks up more from a plant
                
                if (norm_animal_count > 0)
                    % god bless ya OK HERE WE GO
                    % We are looping through all of the animal neighbors,
                    % if we have an animal neighbor, then we are getting
                    % the position of that bad boy.
                    total_population = 0;
                    % If we have an animal, then use the von neighborhood
                    % to get the neighborhood of that cell to determine if
                    % this is the highest pollen cell.
                    for val = 1:length(neighbors)
                        supreme = 0;
                        if (neighbors(val) == ANIMAL)
                            x = von_neighbor_coords(val,1);
                            y = von_neighbor_coords(val,2);
                            % Then loop through all of the neighbors of
                            % this animal, looking for pollen cells, and 
                            % test which is highest value of pollen 

                            animal_neighbors = [extended_grid(x-1, y-1), extended_grid(x-1, y), extended_grid(x-1, y+1),...
                                                extended_grid(x, y-1), extended_grid(x, y), extended_grid(x+1, y+1), ...
                                                extended_grid(x+1, y-1), extended_grid(x+1, y), extended_grid(x+1, y+1)];
                            animal_n_pollen = [extended_pollen_grid(x-1, y-1), extended_pollen_grid(x-1, y), extended_pollen_grid(x-1, y+1),...
                                                extended_pollen_grid(x, y-1), extended_pollen_grid(x, y), extended_pollen_grid(x+1, y+1), ...
                                                extended_pollen_grid(x+1, y-1), extended_pollen_grid(x+1, y), extended_pollen_grid(x+1, y+1)];
                            for n = 1:length(animal_neighbors)
                                if (animal_neighbors(n) == POLLEN)
                                    conc = animal_n_pollen(n);
                                    % Which out of all of these will be the
                                    % right num?
                                    if (conc > supreme)
                                        supreme = von_neighbors(n);
                                    end
                                end
                            end
                            
                            if (current_pollen_conc > supreme)
                                % So if the current pollen cell is the
                                % highest cell, it will become that animal
                                % cell and get its animal concentration,
                                % otherwise this was all for literally
                                % nothing
                                updated_cell = ANIMAL;
                                updated_aniaml_pop = supreme;
                            end
                        end
                    end                              
                    
                elseif (plant_count > 0)
                    % Collecting pollen from all the plants nearby
                    updated_cell = POLLEN;
                    %updated_pollen_cell = current_pollen_conc + pollen_diffuse_rate * plant_n_pollen_conc;                   
                elseif (empty_count > 0)
                    % We are going to lose pollen now, farewell, lovelies
                    % This is loss to all of the surrounding empty cells
                    updated_cell = POLLEN; %Maybe? Maybe not?
                    %updated_pollen_cell = current_pollen_conc - empty_count * pollen_diffuse_rate * current_pollen_conc;                   
                elseif (poll_animal_count > 0)
                    % Pollinated animals are going to snatch this pollen
                    % and take over the entire space, greedy bastards
                    % Wait is this right? Are we sending the pollinated
                    % animal to any pollen cell? IDK now ugh k
                    updated_cell = POLLINATED_ANIMAL;
                    updated_animal_pop_cell = current_animal_pop_cell;
                    %updated_pollen_cell = current_pollen_conc;
                elseif (pollen_count > 0)
                    % We are going to give some pollen to our pollen
                    % friends, bc pollen supports pollen
                    updated_cell = POLLEN;
                    % Actually I don't want to think about this right now
                else 
                    updated_cell = POLLEN;
                    updated_animal_pop_cell = 0;
                    %updated_pollen_cell = current_pollen_conc;
                end
                
            % Ahh so easy that's the good life
            elseif (current_cell == GROWING_PLANT)
                updated_cell = PLANT;
                updated_animal_pop_cell = 0;
                %updated_pollen_cell = 0;
                
            % It's plant time, it just needs a pat pat
            elseif (current_cell == PLANT)
                if (poll_animal_count > 0 && rand < prob_pollination)
                    % So if there are pollinated animals near, they really
                    % want to get rid of their pollen now
                    updated_cell = POLLINATED_PLANT;
                    %updated_pollen_cell = current_pollen_conc + neighbor_poll_animal_conc;
                elseif (rand < prob_plant_death)
                    % Rip a soldier, they tried, but not hard enough
                    updated_cell = EMPTY;
                    %updated_pollen_cell = 0;
                else
                    % These plants are going to produce pollen every dt
                    updated_cell = PLANT;
                    %updated_pollen_cell = current_pollen_conc + init_pollen_count;
                end
                
            elseif (current_cell == POLLINATED_PLANT)
                if (empty_count > 0)
                    % We are going to have to get rid of some pollen sadd
                    % updated_pollen_cell = current_pollen_conc - empty_count * pollen_diffuse_rate * poll_plant_n_conc;
                    updated_cell = POLLINATED_PLANT; % Unless no more pollen!
                elseif (rand < prob_plant_death) 
                    updated_cell = EMPTY;
                    %updated_pollen_cell = 0;
                else
                    % These plants also like to produce pollen every dt
                    updated_cell = POLLINATED_PLANT;
                    %updated_pollen_cell = current_pollen_conc + init_pollen_count;
                end
            end
                
            % Updating next grid with the new cell value
            grids(row-2, col-2 , frame) = updated_cell;
            animal_pop_grids(row-2, col-2, frame) = updated_animal_pop_cell;
            pollen_conc_grids(row-2, col-2, frame) = updated_pollen_cell;
        end
    end
    
    % Recalculate the number of animals and plants for tracking
    animal_counter(frame) = sum(sum(grids(:,:,frame)==ANIMAL)) +...
                              sum(sum(grids(:,:,frame)==POLLINATED_ANIMAL));
    plant_counter(frame) = sum(sum(grids(:,:,frame)==PLANT)) + ...
                             sum(sum(grids(:,:,frame)==GROWING_PLANT)) + ...
                             sum(sum(grids(:,:,frame)==POLLINATED_PLANT));
    animal_pop_counter(frame) = sum(sum(animal_pop_grids(:,:,frame)));
    pollen_conc_counter(frame) = sum(sum(pollen_conc_grids(:,:,frame)));

end

disp("All grids calculated");

%% Visualize the grid
%% Visualize the grid

% Create the window for the animation
viz_fig = figure;
viz_fig2 = figure;

viz_axes = axes(viz_fig);

% Set the colors
map = [ 1       1       1;          % Empty Cell: white
        0  1  0;                    % Growing Plant Cell: Lime Green
        109/255 188/255 0;          % Plant Cell: Green
        109/255 120/255 0;          % Pollinated plant: Camo Green
        125/255 100/255  70/255;    % Animal Cell: Brown 
        0       0      0;           % Pollinated Animal Cell: Black
        230/255 255/255 0];         % Pollen color: Yellow     

colormap(viz_axes, map); 

% Remove axis labels, make aspect ratio look good, and maintain that state
% axis off;
% axis equal;
% hold on;

% viz_fig2 = figure;

disp("Drawing...");
for i = 1:numIterations

    
    % Uncomment following line to allow for frame-by-frame viewing of grid
      w = waitforbuttonpress;
    
    heatmap(viz_fig, grids(:,:,i));
    h=heatmap(viz_fig2, pollen_conc_grids(:,:,i))
%     hmap = heatmap(viz_fig, animal_pop_grids(:,:,i));
    caxis(h, [0 30]);
    % Turn each grid into an image
%     heatmap(viz_fig2, animal_pop_grids(:,:,i));
%     hmap = heatmap(viz_fig, animal_pop_grids(:, :, i));
%     caxis(hmap, [0 30]);
    
    if(i>1) % Refresh the image
        delete(time_counter_text);
        delete(animal_counter_text);
        delete(plant_counter_text);
    end
    
    % Turn each grid into an image
    %image(viz_axes, grids(:, :, i));
    
%     % Draw text
    % time_counter_text = text(viz_axes,0.1,1.025,"Current Time Step: " + ...
    %     string(i), 'Units', 'Normalized');
    % animal_counter_text = text(viz_axes,0.9,0,"Animals: " + ...
    %         animal_counter(i), 'Units', 'Normalized');
    % plant_counter_text = text(viz_axes,0.9,0.05,"Plants: " + ...
    %     plant_counter(i), 'Units', 'Normalized');
    
    pause(1/animation_fps);
end

% Slows down simulation
% s1 = plot(nan,'s', 'color', map(1,:), 'MarkerFaceColor', map(1,:));
% s2 = plot(nan, 's', 'color', map(2,:), 'MarkerFaceColor', map(2,:));
% s3 = plot(nan, 's', 'color', map(3,:), 'MarkerFaceColor', map(3,:));
% s4 = plot(nan, 's', 'color', map(4,:), 'MarkerFaceColor', map(4,:));
% s5 = plot(nan, 's', 'color', map(5,:), 'MarkerFaceColor', map(5,:));
% s6 = plot(nan, 's', 'color', map(6,:), 'MarkerFaceColor', map(6,:));
% s7 = plot(nan, 's', 'color', map(7,:), 'MarkerFaceColor', map(7,:));
% legend([s1, s2, s3, s4, s5, s6, s7], {'Empty', 'Growing Plant', 'Plant'...
%     'Pollinated Plant', 'Animal', 'Pollinated Animal', 'Pollen'}, ...
%     'Location', 'bestoutside');

% Create the graphs for the counters
animal_fig = figure;
animal_axes = axes(animal_fig);

plant_fig = figure;
plant_axes = axes(plant_fig);

plot(animal_axes, 1:numIterations, pollen_conc_counter, 'color', [1 0.5 0.17]);


% plot(animal_axes, 1:numIterations, animal_counter, 'color', [1 0.5 0.17]);
% plot(plant_axes, 1:numIterations, plant_counter);

% title(animal_axes, "Animal Population Throughout Simulation");
% ylabel(animal_axes, "Number of animals (incl. pollen-carrying)");
% xlabel(animal_axes, "Time step");
% ylim(animal_axes, [0, max(animal_counter)+5]); % Extend the maximum y limit

% title(plant_axes, "Plant Population Throughout Simulation");
% ylabel(plant_axes, "Number of plants (incl. growing and pollinated)");
% xlabel(plant_axes, "Time step");
% ylim(plant_axes, [0, max(plant_counter)+5]); % Extend the maximum y limit


disp("Simulation complete!");