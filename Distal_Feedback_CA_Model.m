% CA model for generating anisotropic scale-dependent synthetic landscapes
%
% This is a model loosely based on the CA model from Scanlon et al. 2007,
% however it incorporates a scale-dependent distal negative feedback that
% produces regular patterning (not scale-free).
%
% The script displays a video of model progression and saves the synthetic
% landscape and outputs the image as a .mat file
%
% Scanlon, Todd M., et al. "Positive feedbacks promote power-law clustering 
% of Kalahari vegetation." Nature 449.7159 (2007): 209-212.

clear
clc

format long

%Grid Size
nx=100;

iterations=1000;

%Fractional Tree Cover
ft=.5;

%Choose which model (1=pareto, 2=exponential)
model=1;

%Choose whether to include a distal negative feedback (no = 0 , yes = 1)
dis_feed=1;

%Pareto weighting exponent
p_exp=2;

%Exponential k
k=.2;

%Distal Negative Feedback Range (meters) (must be multiple of cell size)
%This also functions as the window size if dis_feed=0
neg_range=200;

%Pixel length (meters)
cell_size=20;

%Anisotropic Facilitation (for no effect, set upstream=1)
upstream=10;

r=rand(nx,nx);
m=round(r);

z=100;
in=(0:1:z);
[x,y]=meshgrid(in,in);
distance=sqrt(x.^2+y.^2);
distance=distance*cell_size;
sort_distance=sort(distance(:));
pt=zeros(length(sort_distance),1);

if model==1
distfunc=@(x) ((1./(x)).^p_exp);
end

if model==2
distfunc=@(x) exp(-k*x);
end

sort_distance(1)=[];

for j=1:length(sort_distance)
pt(j)=distfunc(sort_distance(j));
end

cum_pt=zeros(length(pt),1);
for j=1:length(cum_pt)
cum_pt(j) = sum(pt(1:j))./sum(pt);
end

cumulative_P = cum_pt(find(sort_distance==neg_range,1,'first'))

dist_window=distance(1:neg_range/cell_size+1,1:neg_range/cell_size+1);
dist_window(dist_window>neg_range)=0;

neg_window=dist_window;
neg_window(neg_window>0)=1;
neg_window(1)=1;
neg_window=[fliplr(neg_window) neg_window(:,2:end)];
neg_window=[flipud(neg_window);neg_window(2:end,:)];

window_size=length(neg_window);
neighb_edge=(window_size+1)/2;
neighb_edge2=neighb_edge-1;

a=distfunc(dist_window);
a(a==inf | a==1)=0;
a=[fliplr(a) a(:,2:end)];
a=[flipud(a);a(2:end,:)];
suma=sum(a(:));
a(1:end,neighb_edge)=a(1:end,neighb_edge)*upstream;
a=a.*suma/sum(a(:));

i_matrix=ones(size(m));
i_matrix=padarray(i_matrix,[neighb_edge2 neighb_edge2]);
i_column=im2col(i_matrix,[window_size window_size],'sliding');

a_column=reshape(a,numel(a),1);
a_column=repmat(a_column,1,nx*nx);
a_column=a_column.*i_column;

neg_column=reshape(neg_window,numel(neg_window),1);
neg_column=repmat(neg_column,1,nx*nx);
neg_column=neg_column.*i_column;
neg_total=sum(neg_column);

pt_fun=@(m_column,i_column) sum(m_column.*a_column)./sum(a_column);

count=0;
for n=1:iterations
m_pad=padarray(m,[neighb_edge2 neighb_edge2]);
m_pad(:,1)=1;
m_pad(:,end)=1;
m_column=im2col(m_pad,[window_size window_size],'sliding');

pt_column=pt_fun(m_column,i_column);
pt_matrix=reshape(pt_column,nx,nx);

if dis_feed == 1
    ft_local_column=sum(m_column.*neg_column)./neg_total;
    ft_local=reshape(ft_local_column,nx,nx);
    ft_actual=ft_local;
else
    ft_actual=sum(sum(m))/numel(m);
end

Pot=pt_matrix+(ft-ft_actual)./(1-ft_actual);
Pto=1-pt_matrix+(ft_actual-ft)./ft_actual;
rand1=rand(size(m));
rand2=rand(size(m));
trans_ot=(.2>=rand1 & m(:,:)==0 & Pot>=rand2);
%Pto=1-sum(sum(trans_ot.*pt_matrix))/sum(sum(trans_ot))+(ft_actual-ft)/ft_actual;
trans_to=(.2>=rand1 & m(:,:)==1 & Pto>=rand2(:,:));

m=m+trans_ot;
m=m-trans_to;

count=count+1;
if count==1000
    disp(n)
    count=0;
end

% Shows a video of the model progression
 imshow(1-m);
 drawnow

end;

imshow(1-m);

save('serpentine.mat','m')