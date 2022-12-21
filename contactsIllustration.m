% script contacts.m
% author: T. Nguyen-Huu 
%
% description: this script generates 'iter' iteration of random populations
% of 'N' individuals, and computes the number of infection when there are
% 'I' in 'Irange' infected individuals, for an average contact number of 
% 'beta' in 'betaRange' per individual, (1) without isolation and (2) for a 
% proportion 0<'v'<1 of isolated people.
% The values obtained are compared with the theoretical values beta.IS/N
% and beta.(1-v)^2.IS/N for each run.


N = 20;
Irange = [floor(N/10), 7*floor(N/10)];
betaRange = [2 4];


% generation of 'iter' cases with random values for 'v', 'beta' and initial
% condition 'I'
v = randi(N);
beta = betaRange(1)+rand*(betaRange(2)-betaRange(1));
I = randi(Irange);
[before, after]= count_contacts(N,I,beta,v);


% The function 'count_contacts' generate a random population of 'N' individuals
% with 'I' infected individuals, and create a contact graph such that the
% average number of contacts for each individual is 'beta' (we considered 
% that the probability of infection p=1 for the sake of simplicity). 
% It returns the number of infections occuring because of the contacts 
% before isolation ('nb_infections'), and after setting a containment of 
% 'nb_removed' individuals ('nb_infections_after_containment').
% The values obtained are compared with the theoretical values beta.IS/N
% and beta.(1-v)^2.IS/N.

function [nb_infections, nb_infections_after_containment] = count_contacts(N,I,beta,nb_removed)
    blue = [41, 128, 185]/255;
    red = [192, 57, 43]/255;
    grey = [100, 100, 100]/255;
    light_grey = [200, 200, 200]/255;

    % create a contact list with N*beta/2 elements
    contact_list = cell(floor(N*beta/2),1);
    [contact_list{1:end}] = deal([0,0]);
    
    % populate the contact_list with links [i,j] such as i ~= j. For
    % practical aspect, links are arranged such that i<j.
    for i=1:numel(contact_list)
        contact_added = false;
        while ~contact_added 
            person1 = randi(N);
            person2 = person1;
            while person2 == person1
                person2 = randi(N);
            end
            contact_edge = [min(person1,person2), max(person1,person2)];
            contact_added = ~iscellmember(contact_edge,contact_list);
            if contact_added
                contact_list{i} =contact_edge;
            end
        end
    end

    % chose I infected individuals
    I_list = randperm(N,I);

    % computation of the number of infectious links (without isolation)
    nb_infections = sum(contact_is_infectious(contact_list,I_list));
    fprintf(['No containment: beta.IS/N = %4.2f. Number of infectious contacts ' ...
        ': %d (total contacts: %d).\n'], beta*I*(N-I)/N, nb_infections,numel(contact_list));
    
    % put 'nb_removed' individuals in isolation and remove the contact links involving isolated individuals
    to_be_removed = cellfun(@(x)(x(1)<=nb_removed),contact_list); 
    contact_list_after_containment = contact_list(~to_be_removed,:);

     % computation of the number of infectious links (after isolation)
    nb_infections_after_containment=sum(contact_is_infectious(contact_list_after_containment,I_list));
    fprintf(['With containment: beta.(1-v)².IS/N = %4.2f. Number of infectious ' ...
        'contacts: %d (%d individuals were isolated).\n\n'], ...
        (1-nb_removed/N)^2*beta*I*(N-I)/N, nb_infections_after_containment,nb_removed);

    s = cellfun(@(x)x(1),contact_list);
    t = cellfun(@(x)x(2),contact_list);

    infectious_edges = contact_is_infectious(contact_list,I_list);
    G = graph(s,t,table(infectious_edges,'VariableNames',{'Infected'}));
    
    color_list = repmat(blue,N,1);
    color_list(I_list,:) = repmat(red,I,1);

    edge_color_list = repmat(grey,numel(contact_list),1);
    edge_color_list(G.Edges.Infected,:) = repmat(red,nb_infections,1);
    
    style_list = repmat((":"),numel(contact_list),1);
    style_list(G.Edges.Infected,:) = repmat(("-"),nb_infections,1);

    inactivated_edges = cellfun(@(x)(x(1) <= nb_removed),contact_list); 
    G2 = graph(s,t,table(infectious_edges, inactivated_edges,'VariableNames',{'Infected','Inactivated'}));
    to_deactivate = and(G2.Edges.Infected,G2.Edges.Inactivated);

    style_list_after = style_list;
    edge_color_list_after = edge_color_list;
    edge_color_list_after(to_deactivate,:) = repmat(grey,sum(to_deactivate),1);

    style_list_after(to_deactivate,:) = repmat(("--"),sum(to_deactivate),1);

    label_colors = repmat([0 0 0],N,1);
    label_colors(1:nb_removed,:) = repmat(light_grey,nb_removed,1);

    fig = figure(1);
    WINDOW_SIZE = [1200 600];
    Pix_SS = get(0,'screensize');
    fig.Position = [(Pix_SS(3)-WINDOW_SIZE(1))/2   (Pix_SS(4)-WINDOW_SIZE(2))/2   WINDOW_SIZE(1) WINDOW_SIZE(2)];
    fig.Color = [1 1 1];
    %     clf;
%     hold on;
    fig1 = subplot(1,2,1);
    plot(G,'Layout','circle','NodeColor',color_list,'EdgeColor',edge_color_list, ...
        'LineStyle',style_list);
    fig1.OuterPosition = [0.01 0.01 0.5 1];
    fig1.Position = [0.01 0.01 0.5 1];
    fig1.Box = 'off';
    axis off;
    fig2 = subplot(1,2,2);
    plot(G2,'Layout','circle','NodeColor',color_list,'EdgeColor', ...
        edge_color_list_after,'LineStyle',style_list_after,'NodeLabelColor',label_colors);
    fig2.OuterPosition = [0.5 0 0.5 1];
    fig2.Position = [0.5 0 0.5 1];
    fig2.Box = 'off';
    axis off;
end

% Auxiliary function similar to 'ismember', but for cell arrays instead of
% arrays. Returns 1 if one of the cells of 'cellArray' is equal to 'elem', 
% 0 otherwise.

function res = iscellmember(elem,cellArray)
    res = any(cellfun(@(x)all(x == elem),cellArray));
end

% Returns a cell array 'res', such that the kth cell is 1 is the kth 
% contact in 'ContactsCellArray' is infectious, and 0 otherwise. Each cell 
% of 'ContactsCellArray' contains an contact link [i,j] between two 
% individuals i and j. If i is Infected (i is in 'I_list') and j Susceptible 
% (j is not in 'I_list',or inversely, then the contact is infectious. 

function res = contact_is_infectious(contactsCellArray,I_list)
    res = cellfun(@(x)xor(ismember(x(1),I_list),ismember(x(2),I_list)),contactsCellArray);
end
