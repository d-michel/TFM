%%%%%%%%%%  DENSITY  %%%%%%%%%%
data = readtable("density_16.txt");
energy = data.(1);
nstates = data.(2);

hold on
set(gca, 'fontsize', 13)
box on
figure(1)
plot(energy./16, nstates./sum(nstates), "-ob", linewidth=1.5)
xline(1.1,":g",linewidth=1.5)
xlabel("$E/N$", "Interpreter", "latex",'fontsize', 18)
ylabel("$D(E)$", "Interpreter", "latex",'fontsize', 18)

%%%%%%%%%%  IPR  %%%%%%%%%%
data = readtable("IPR_150_10_1_6.txt");
energy = data.(1);
IPR = data.(2);

hold on
set(gca, 'fontsize', 13)
box on
figure(2)
plot(energy./10, IPR, ".b", linewidth=1.5)
xline(1.1,":g",linewidth=1.5)
xlabel("$E/N$", "Interpreter", "latex",'fontsize', 18)
ylabel("IPR", "Interpreter", "latex",'fontsize', 18)

%%%%%%%%%%  PERES LATTICES %%%%%%%%%%
data = readtable("peres_lattices_150_10_1_6_sign.txt");
energy = data.(1);
s00 = data.(2);
s11 = data.(3);
s22 = data.(4);
ada = data.(5);
a_plus_ad = data.(6);

hold on
set(gca, 'fontsize', 13)
box on
figure(3)
scatter(energy./10,s00,50,a_plus_ad,".")
colorbar
xline(1.1,":g",linewidth=1.5)
xlabel("$E/N$","Interpreter","latex",'fontsize', 18)
ylabel("$\langle\hat{\Sigma}_{00}\rangle$","Interpreter","latex",'fontsize', 18)

hold on
set(gca, 'fontsize', 13)
box on
figure(4)
scatter(energy./10,s11,50,a_plus_ad,".")
colorbar
xline(1.1,":g",linewidth=1.5)
xlabel("$E/N$","Interpreter","latex",'fontsize', 18)
ylabel("$\langle\hat{\Sigma}_{11}\rangle$","Interpreter","latex",'fontsize', 18)

hold on
set(gca, 'fontsize', 13)
box on
figure(5)
scatter(energy./10,s22,50,a_plus_ad,".")
colorbar
xline(1.1,":g",linewidth=1.5)
xlabel("$E/N$","Interpreter","latex",'fontsize', 18)
ylabel("$\langle\hat{\Sigma}_{22}\rangle$","Interpreter","latex",'fontsize', 18)

hold on
set(gca, 'fontsize', 13)
box on
figure(6)
scatter(energy./10,ada,50,a_plus_ad,".")
colorbar
xline(1.1,":g",linewidth=1.5)
xlabel("$E/N$","Interpreter","latex",'fontsize', 18)
ylabel("$\langle\hat{a}^\dagger\hat{a}\rangle$","Interpreter","latex",'fontsize', 18)

%%%%%%%%%%  STATISTICS  %%%%%%%%%%
% data = readtable("statistics_150_10_1_6.txt");
% energy = data.(1);
% ratio = data.(2);
% min_ratio = data.(3);
% figure(1)
% hold on
% xlim([0 50])
% histogram(ratio,Normalization="pdf")
% Pr_GOE = @(r) (96/25) * (r + r^2)/((1 + r)^2 - 4*r/5)^(5/2);
% Pr_Poisson = @(r) 1/(1 + r)^2;
% fplot(Pr_GOE,[0 50])
% fplot(Pr_Poisson,[0 50])
% figure(2)
% hold on
% histogram(min_ratio,Normalization="pdf")
% Pr_GOE = @(r) 2 * (96/25) * (r + r^2)/((1 + r)^2 - 4*r/5)^(5/2);
% Pr_Poisson = @(r) 2/(1 + r)^2;
% fplot(Pr_GOE,[0 1])
% fplot(Pr_Poisson,[0 1])
% mean_window_min_ratio = fast_window_analysis(energy,min_ratio,15);
% x = (mean_window_min_ratio(:,2) + mean_window_min_ratio(:,3))./20;
% figure(3)
% plot(x,mean_window_min_ratio(:,1),"o")
% mean(ratio)
% mean(min_ratio)
% 
% function mean_window_min_ratio = fast_window_analysis(energy, min_ratio, num_windows)
%     end_points = ones(num_windows+1,1);
%     length_window = (energy(length(energy)) - energy(1))/num_windows;
%     mean_window_min_ratio = zeros(num_windows,3);
%     j = 1;
%     for i = 1:num_windows
%         while (energy(j) < (energy(1) + (i+1)*length_window)) && (j < length(energy))
%             j = j + 1;
%         end
%         end_points(i+1) = j;
%         mean_window_min_ratio(i,1) = mean(min_ratio(end_points(i):end_points(i+1)));
%         mean_window_min_ratio(i,2) = energy(end_points(i));
%         mean_window_min_ratio(i,3) = energy(end_points(i+1));
%     end
% end

%%%%%%%%%%  SYMMETRIES  %%%%%%%%%%
str = ["symmetries_150_10_1_6" "50_10_1_12_f1_E10-4"];
for i = 1:length(str)
    data = readtable(str(i)+".txt");
    energy = data.(1);
    operator = data.(2);
    index = find(energy./10>6,1);
    %index = length(operator);

    hold on
    set(gca, 'fontsize', 13)
    box on
    figure(i+6)
    plot(energy(1:index)./10, operator(1:index), ".b", linewidth=1.5)
    xlim([-2 7])
    xline(1.1,":g",linewidth=1.5)
    xlabel("$E/N$", "Interpreter", "latex",'fontsize', 18)
    ylabel("$\langle \hat{a}+\hat{a}^\dagger\rangle$", "Interpreter", "latex",'fontsize', 18)
end

%%%%%%%%%%  Ec  %%%%%%%%%%
data = readtable("Ec_150_spectrum_10.txt");
energy = data.(1);
Ec = data.(2);

hold on
set(gca, 'fontsize', 13)
box on
figure(10)
plot(energy./i, Ec, ".b", linewidth=1.5)
xline(1.1,":g",linewidth=1.5)
xlabel("$E/N$", "Interpreter", "latex",'fontsize', 18)
ylabel("$\langle\hat{\mathcal{C}}\rangle$", "Interpreter", "latex",'fontsize', 18)
hold on
set(gca, 'fontsize', 13)
box on