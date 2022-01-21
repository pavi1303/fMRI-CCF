cd('D:\LRCBH\Final')
addpath(genpath('C:\Users\pavig\Documents\CCF Project\fMRI-multimodal\Multivariate decomposition\ICA'))
pca_tcat1 = uint8(pca_tcat1);
pca_tcat2 = uint8(pca_tcat2);
pca_tcat3 = uint8(pca_tcat3);
pca_tcat4 = uint8(pca_tcat4);
pca_tcat5 = uint8(pca_tcat5);
pca_tcat6 = uint8(pca_tcat6);
pca_tcat7 = uint8(pca_tcat7);
pca_tcat8 = uint8(pca_tcat8);
pca_tcat9 = uint8(pca_tcat9);
pca_tcat10 = uint8(pca_tcat10);
pca_tcat11 = uint8(pca_tcat11);

pca_tcat_1 = vertcat(pca_tcat1,pca_tcat2,pca_tcat3,pca_tcat4,pca_tcat5);
pca_tcat_2 = vertcat(pca_tcat6,pca_tcat7,pca_tcat8,pca_tcat9,pca_tcat10,pca_tcat11);
pca_tcat = vertcat(pca_tcat_1, pca_tcat_2);
pca_tcat = uint8(pca_tcat);

[S1,W1,White1,E1,eigval1,convergence1,A1,B1,A_reduced1,X_reduced1,Sigma_reduced1]=ica_DC_improved(pca_tcat,0,11,1E-6,30,[],1,1,[]);
