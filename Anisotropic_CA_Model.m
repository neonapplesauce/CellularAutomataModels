% Cellular automata model for generating anisotropic synthetic landscapes
%
% This model is loosely based on the CA model from Scanlon et al. 2007,
% however expands the model to allow for an anisotropic kernel (one where
% the weighting of the positive feedback is greater in one direction).
% Similar to the model proposed by Acharya et al., this model was originally
% developed to simulate the Everglades ridge slough landscape.
%
% Importantly, this model can also simulate two degradation mechanisms: one 
% with subtractive degradation (patches are randomly destroyed, such as may 
% be the case with rapid increased hydroperiod) and one with additive
% degradation (such as may be the case with rapid decreased hydroperiod).
%
% The script also produces a video of landscape progression and outputs
% cumulative distribution functions (CDF) of patch size, as well as a plot
% of patch perimiter vs area (slope of plot is fractal dimension).
%
% References:
% Acharya, Subodh, et al. "Coupled local facilitation and global hydrologic 
% inhibition drive landscape geometry in a patterned peatland." Hydrology and
% Earth System Sciences 19.5 (2015): 2133-2144.
%
% Scanlon, Todd M., et al. "Positive feedbacks promote power-law clustering 
% of Kalahari vegetation." Nature 449.7159 (2007): 209-212.

clear
clc

tstart=tic;

%Grid Size
nx=200;

%pixelsize   %Changing this parameter doesn't seem to do anything
dmin=1;

% Which kernel model to use?  1=pareto 2=exponental 3=linear
model=1;

% Subtractive degradation level. Ranges from 0-1 (0=none, 1=max) 
s_degradation = .1;

% Additive degradation level. Ranges from 0-1 (0=none, 1=max)
a_degradation = .1;

iterations = 10000;

%Fractional Tree Cover
ft=.5;

%Elongation (1=isotropic)
elongation=1.5;

%Cumulative P accounted for by the distance weighting neighborhood
neighb_fraction=.9999;

%Pareto weighting exponent
k=7;

%Proportion of cells to transition      %seems to change curvature in lowerend of distrubtion 
tamount=.1;                             %(high values mean lots of small patches)
                             
%random 50/50 initiation
r=rand(nx,nx);
m=logical(round(r));

z=100;
in=(0:1:z);
[x,y]=meshgrid(in,in);
y=y/elongation;
distance=sqrt(x.^2+y.^2);
sort_distance=sort(distance(:));
pt=zeros(length(sort_distance),1);

switch model
    case 1
        for j=1:length(sort_distance)
        pt(j)=((dmin./(sort_distance(j))).^k);
        end;
    case 2
        for j=1:length(sort_distance)
        pt(j)=exp(-k*sort_distance(j));
        end;
    case 3
        for j=1:length(sort_distance)
        pt(j)=(1-sort_distance(j)/k);
        if pt(j)<0
            pt(j)=0;
        end
        end;
end

cum_pt = 1-(pt(2:end)./sum(pt(2:end)));
cum_pt_thresh=cum_pt(cum_pt>=neighb_fraction);

thresh_ind=length(cum_pt)-length(cum_pt_thresh);
dis=sort_distance(thresh_ind+1);

% neighb_edge2=neighb_edge-1;
distance0=distance;
distance0(distance0>dis)=0;
kernel=distance0(1:find(max(distance0')==0,1,'first')-1,1:find(max(distance0)==0,1,'first')-1);
kernel=[fliplr(kernel) kernel(:,2:end)];
kernel=[flipud(kernel);kernel(2:end,:)];
window_size=size(kernel);
neighb_edge=(window_size-1)/2;

switch model
    case 1     
        a=(dmin./(kernel)).^k;
    case 2
        a=exp(-k.*kernel);
    case 3
        a=(1-kernel./k);
end

a(a==inf)=0;
a_all=a(:);

i_matrix=ones(size(m));
i_matrix=padarray(i_matrix,neighb_edge);
i_column=im2col(i_matrix,window_size,'sliding');
i_column=logical(i_column);

z=(i_column'*a_all)';

clearvars -except a_all z iterations m nx ft tamount neighb_edge neighb_edge2 window_size a_degradation s_degradation

% pt_fun=@(m_column) sum(m_column.*a_column)./z;
pt_fun=@(m_column) (m_column'*a_all)'./z;

count=1;
for n=1:iterations
    % disp(n)

m_pad=logical(padarray(m,neighb_edge));
m_column=im2col(m_pad,window_size,'sliding');

pt_column=pt_fun(m_column);
pt_matrix=reshape(pt_column,nx,nx);

ft_actual=sum(sum(m))/numel(m);

Pot=pt_matrix+(ft-ft_actual)/(1-ft_actual);
Pto=1-pt_matrix+(ft_actual-ft)/ft_actual;

rand1=rand(size(m));
rand2=rand(size(m));
trans_ot=(tamount>=rand1 & m(:,:)==0 & Pot>=rand2);
trans_to=(tamount>=rand1 & m(:,:)==1 & Pto>=rand2);

m=m+trans_ot;
m=m-trans_to;
m=logical(m);

% Shows a video of the landscape progression
if count==20
figure(1)
    imshow(1-m);
    drawnow
    count=0;
    
% Shows a video of patch size CDF progression
% c=bwlabel(m,4);
% rstats=regionprops(c,'all');
% areas=cat(1,rstats.Area);
% figure(2);
% loglog(unique(areas),ccdfcalc(areas),'.');
% drawnow

end
count=count+1;
end;
title('Original Landscape')

% Simulate degradation (adds in random noise)
msave=m;
% subtractive degradation
m=m-(rand(size(m))>(1-a_degradation));
% additive degradation
m=m+(rand(size(m))>(1-s_degradation));
m=m>0;
figure(2)
imshow(1-m); 
title('Degraded Landscape')

% Plot final CDF of patch areas  
c=bwlabel(m,4);
rstats=regionprops(c,'all');
areas=cat(1,rstats.Area);
figure(3)
loglog(unique(areas),ccdfcalc(areas),'.')
title('Patch size CDF (degraded)')

% Plot perimeter vs patch area
        data_pad=padarray(m,[1 1]);    
        cellperim=zeros(size(data_pad));
        for j=2:nx+1 
            for k=2:nx+1
                if data_pad(j,k)==1
                    cellperim(j,k)=sum([data_pad(j-1,k) data_pad(j+1,k) data_pad(j,k+1) data_pad(j,k-1)]==0);
                end
            end
        end
        cellperim=cellperim(2:nx+1,2:nx+1);
        perims=sum(sum(cellperim));
        
        perimindex_p=cellperim(cellperim~=0);
        perimindex_c= c(cellperim~=0);
        
        ridgeperims=zeros(length(areas),1);
        for n1=1:length(areas)
            ridgeperims(n1)=sum(perimindex_p(perimindex_c==n1));
        end
        figure(4) 
        loglog(ridgeperims,areas,'.')
        title('Perimiter vs Area (degraded)')
        
% Save final landscape image        
save('anisotropic_landscape.mat','m')