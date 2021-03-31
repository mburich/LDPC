clear all
close all
clc



%% BSC plots


figure
pflip = [0.05       0.05368421 0.05736842 0.06105263 0.06473684 0.06842105 0.07210526 0.07578947 0.07947368 0.08315789 0.08684211 0.09052632 0.09421053 0.09789474 0.10157895 0.10526316 0.10894737 0.11263158 0.11631579 0.12      ];
ber100pckts = [0.        0.        0.        0.        0.        0.        0. 0.        0.0014775 0.023816  0.0534875 0.0680405 0.0756385 0.083936 0.0901925 0.0966725 0.1021875 0.107345  0.1125035 0.116881 ];
deber = [0.         0.         0.         0.         0.         0.  0.         0.         0.         0.         0.0570185  0.06753978  0.07589482 0.08315725 0.08970529 0.09573177 0.10134833 0.10663002  0.11164337 0.11640694];
capacity = 0.11002;
semilogy(pflip,ber100pckts,'--rx'); hold on
semilogy(pflip,max(eps,deber),'--ko');

xlabel('Channel flip probability');
ylabel('BER');
title(['BSC - Regular (3,6) - Capacity at ',num2str(capacity)])
grid on
legend(['BER - 100 packets - N=10000'],['Density Evolution'],'location','SouthEast')
ylim([1e-4, 4e-1])
saveas(gcf,'bsc.png')
close all
%% AWGN plots


figure

esn0 = [0.8        0.83157895 0.86315789 0.89473684 0.92631579 0.95789474 0.98947368 1.02105263 1.05263158 1.08421053 1.11578947 1.14736842 1.17894737 1.21052632 1.24210526 1.27368421 1.30526316 1.33684211 1.36842105 1.4       ];
ber100pckts = [0.0837705 0.081332  0.079241  0.076574  0.0731045 0.0698125 0.067885 0.0570835 0.0483585 0.0435565 0.027729  0.0225205 0.013895  0.007919 0.007191  0.0015815 0.000496  0.        0.        0.       ];
deber = [0.08433212 0.08201562 0.0796011  0.077069   0.0743919  0.07152926 0.06841645 0.06493827 0.06084632 0.05531958 0.         0. 0.         0.         0.         0.         0.         0. 0.         0.        ];
capacity = 0.0;
semilogy(esn0,ber100pckts,'--rx'); hold on
semilogy(esn0,max(eps,deber),'--ko');

xlabel('EbN0 [dB]');
ylabel('BER');
title(['AWGN - Regular (3,6) - Capacity at ',num2str(capacity)])
grid on
legend(['AWGN - 100 packets - N=10000'],['Density Evolution'],'location','SouthEast')
ylim([1e-4, 4e-1])
saveas(gcf,'awgn.png')
close all
