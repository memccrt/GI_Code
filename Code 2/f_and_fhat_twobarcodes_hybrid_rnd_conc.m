function [fhatmat,fmat]=f_and_fhat_twobarcodes_hybrid_rnd_conc(probes,R,noise,ntrials)
%inputs: probes - fluorescent proteins combinations that make up each probe 
%        R - emission spectra of each probe for each excitation channel 
%        noise - standard deviation of normal distribution (zero mean)
%        ntrials - number of replicates for each barcode
%output: fhatmat - predicted propbes present in mixture
%        fmat - actual probes present in mixture
%Purpose: generate f and fhat data for each barcode

N_PROBES=length(probes(:,1));
N_BATCH_A=round(N_PROBES/2);
n=ntrials;
lower_conc = 0.5;
higher_conc = 1.5;
for i=1:n
    parfor j = 1:N_PROBES
    rng('shuffle');
        f = zeros(N_PROBES,1);
        f(j)=1;
        if j<= N_BATCH_A
            f(randi([1 N_BATCH_A],1,1))=lower_conc+(higher_conc-lower_conc)*rand;
            f(randi([N_BATCH_A+1 N_PROBES],1,2))=lower_conc+(higher_conc-lower_conc)*rand;
        else
            f(randi([N_BATCH_A+1 N_PROBES],1,1))=lower_conc+(higher_conc-lower_conc)*rand;
            f(randi([1 N_BATCH_A],1,2))=lower_conc+(higher_conc-lower_conc)*rand;
        end
        mu = R * f ; %calculate mu from ref matrix and f
        Length_mu = length(mu); 
        temp=noise.*randn(Length_mu,1); %creates variable for noise
        mu = mu + temp; %add noise to mu
        f_hat = Unmix(R,mu,0,zeros(N_PROBES,1),[],1);
        fmat(:,j,i) = f;
        fhatmat(:,j,i) = f_hat;
        f_hat = zeros(N_PROBES,1);
    end
end
fmat=reshape(fmat,N_PROBES,N_PROBES*n);
fhatmat=reshape(fhatmat,N_PROBES,N_PROBES*n);