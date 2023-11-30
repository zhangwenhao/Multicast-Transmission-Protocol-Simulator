clc;
clear;
close all;

parpool(10);

% Radius of the area
radius = 1000;

% Target rate
R = 0.1;
% Channel usage count
N = 100;
% Data size
M = 1 * 1024 * 1024;

% Noise variance
sigma_k = fun_db_to_math(-70);
sigma_w = sigma_k;
% UAV height
h = 500;
% Channel coefficient
beta = fun_db_to_math(-40);
% Transmission probability
p1 = 0.5; 
% Covertness constraint
epsilon = 0.1;


looptimes=20;
POINT_TIME=zeros(1, looptimes, 'double');
MG_TIME=zeros(1, looptimes, 'double');
UNIFORM_TIME=zeros(1, looptimes, 'double');

loop_times = zeros(1,looptimes,'double');

varsigma = sqrt(N/(2*pi*(exp(2*R)-1)));
vartheta = exp(R)-1;
eta_k = @(gamma_ak_temp) -varsigma*(gamma_ak_temp-vartheta)+0.5;
C_k = @(eta_k_temp) N*R*(1-eta_k_temp);

warden=[0,0];
r_w = 150;

numberNodes=200;
cluster_parameter=0.5;

distributionLoop=10000;
groupLoop=20;

parfor index=1:looptimes
    disp(index);
    loop_times(index)=index;
    distribution_Time_p=0;
    distribution_Time_kmeans=0;
    distribution_Time_mbs=0;
    for d_l=1:distributionLoop
        fprintf('distributionloop: %d\n', d_l);
        centerszie_p = poissrnd(index);
        [data, ~]=poisson_cluster(numberNodes, centerszie_p, cluster_parameter, radius); 
        k_p = size(data, 1);
        close all;
        plot_gu(warden,data,radius);
        currentgroup=0;
        for i= 1:k_p
            h_ak=beta/power(h,2);
            h_aw=beta/(power(h,2)+power(norm(data(i,:) - warden),2)); 
            p_ak=4*epsilon*sigma_w*sqrt(2/N)/h_aw;
            gamma_ak = p_ak*h_ak/sigma_k;
            R_eta_k=eta_k(gamma_ak);
            R_C_k=C_k(R_eta_k);
            currentgroup=currentgroup+ M/(p1*R_C_k);
        end
        distribution_Time_p=distribution_Time_p + currentgroup;
        MG_CONS_Time_kmeans=0;
        MG_CONS_Time_mbs=0;
        
        for g_l=1:groupLoop
            fprintf('groupLoop: %d\n', g_l);
            [MBSLocations,finalRadius, sortedWdx]=group_uniform_radius(data,warden,r_w, radius/2);
            %close all;
            %plot_uniform_group(warden, data, MBSLocations, sortedWdx, r_w, radius, finalRadius);

            [idx, ctr, wdx] = group_k_means( data, warden, r_w);
            %plot_mg(warden, wdx, ctr, radius, r_w, idx);

            MG_kmeans = [ctr; wdx];
            MG_mbs = [MBSLocations; sortedWdx];
    
            k_kmeans = size(MG_kmeans, 1);
            k_mbs = size(MG_mbs, 1);
            current_MG_Time_kmeans=0;
            current_MG_Time_mbs=0;
            for i= 1:k_kmeans
                h_ak=beta/(power(h,2)+power(MG_kmeans(i,3),2));
                h_aw=beta/(power(h,2)+power(norm(MG_kmeans(i, 1:2) - warden),2)); 
                p_ak=4*epsilon*sigma_w*sqrt(2/N)/h_aw;
                gamma_ak = p_ak*h_ak/sigma_k;
                R_eta_k=eta_k(gamma_ak);
                R_C_k=C_k(R_eta_k);
                current_MG_Time_kmeans=current_MG_Time_kmeans + M/(p1*R_C_k);
            end
            for i = 1:k_mbs
                h_ak=beta/(power(h,2)+power(MG_mbs(i, 3),2));
                h_aw=beta/(power(h,2)+power(norm(MG_mbs(i, 1:2) - warden),2)); 
                p_ak=4*epsilon*sigma_w*sqrt(2/N)/h_aw;
                gamma_ak = p_ak*h_ak/sigma_k;
                R_eta_k=eta_k(gamma_ak);
                R_C_k=C_k(R_eta_k);
                current_MG_Time_mbs=current_MG_Time_mbs + M/(p1*R_C_k);
            end
            MG_CONS_Time_kmeans=MG_CONS_Time_kmeans+current_MG_Time_kmeans;
            MG_CONS_Time_mbs=MG_CONS_Time_mbs+current_MG_Time_mbs;
        end
        distribution_Time_kmeans=distribution_Time_kmeans+MG_CONS_Time_kmeans/groupLoop;
        distribution_Time_mbs=distribution_Time_mbs+MG_CONS_Time_mbs/groupLoop;
    end 
     POINT_TIME(index) = distribution_Time_p/distributionLoop;
     MG_TIME(index) = distribution_Time_kmeans/distributionLoop;
     UNIFORM_TIME(index) = distribution_Time_mbs/distributionLoop;
end

delete(gcp);

figure(1); 
plot(loop_times,POINT_TIME,'Color','k','LineStyle','--','Marker','*','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:20:length(loop_times),'linewidth',1);
hold on;
plot(loop_times,UNIFORM_TIME,'Color','k','LineStyle',':','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:20:length(loop_times),'linewidth',1);
plot(loop_times,MG_TIME,'Color','k','LineStyle','-','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:20:length(loop_times),'linewidth',1);

xlabel('Mean number of clusters, $\lambda$','Interpreter','latex');
ylabel('Transmission time (ms)');
xlim([1 20]);
xticks(1:1:20);
legend('Point-to-Point','MGs with SMBSP','MGs with K-means++','Location','best');
box on;
grid on;

