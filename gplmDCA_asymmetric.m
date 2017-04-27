% Copyright 2014 - Christoph Feinauer and Marcin J. Skwark (christophfeinauer@gmail.com, marcin@skwark.pl)
% Copyright 2012 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
% All rights reserved
% 
% Permission is granted for anyone to copy, use, or modify this
% software for any uncommercial purposes, provided this copyright 
% notice is retained, and note is made of any changes that have 
% been made. This software is distributed without any warranty, 
% express or implied. In no event shall the author or contributors be 
% liable for any damage arising out of the use of this software.
% 
% The publication of research using this software, modified or not, must include an 
% appropriate citation to:
%       C. Feinauer, M.J. Skwark, A. Pagnani, E. Aurell, Improving Contact Prediction
%       along Three Dimensions. PLoS Comput Biol 10(10): e1003847.
%
% 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
% 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013) 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function gplmDCA_asymmetric(fastafile,outputfile, lambda_h,lambda_J,lambda_G,reweighting_threshold,nr_of_cores,M)
    options.method='lbfgs';	%Minimization scheme. Default: 'lbfgs', 'cg' for conjugate gradient (use 'cg' if out of RAM).
    options.optTol=1e-5; 	%default: 1e-5
    options.progTol=1e-9;   	%default: 1e-9
    options.Display='off';

%    addpath(genpath(pwd))

    
%Read inputfile (removing inserts), remove duplicate sequences, and calculate weights and B_eff.
    [N,B_with_id_seq,q,Y]=return_alignment(fastafile);
    Y=unique(Y,'rows');
    [lH,rH]=gapMat(int32(Y-1));
    lH=int32(lH);
    rH=int32(rH);
    [B,N]=size(Y);
    weights = ones(B,1);
    if M>N
	error('ERROR: Maximal Gap length must be smaller than length of alignment!');
    end	
    if M==-1
        [~,gapHist]=gapCount(int32(Y-1));
        M=max(find(gapHist~=0))+1;
        if M>N-1
                M=N-1;
        end
    end
  
    nrGapParam=M*(N-(M+1)/2+1); 

    disp(['Maximum gap length: ' num2str(M)]);
  
    if reweighting_threshold>0.0
        fprintf('Starting to calculate weights \n...');
        tic
        Y=int32(Y);
        m=calc_inverse_weights(Y-1,reweighting_threshold);
        weights=1./m;

        fprintf('Finished calculating weights \n');
        toc
    end
    B_eff=sum(weights);
    fprintf('### N = %d B_with_id_seq = %d B = %d B_eff = %.2f q = %d\n',N,B_with_id_seq,B,B_eff,q);
   	
%Prepare inputs to optimizer.
    field_lambda=lambda_h*B_eff;   
    coupling_lambda=lambda_J*B_eff/2;    %Divide by 2 to keep the size of the coupling regularizaion equivalent to symmetric variant of plmDCA.
    gap_lambda=lambda_G*B_eff/N; %Regularization is added for every g_r!
    Y=int32(Y);q=int32(q);
    w=zeros(q+q^2*(N-1)+nrGapParam,N); %Matrix in which to store parameter estimates (column r will contain estimates from g_r).
%Run optimizer.
    if nr_of_cores>1
        matlabpool('open',nr_of_cores)   
        tic
        parfor r=1:N
            disp(strcat('Minimizing g_r for node r=',int2str(r)))       
            wr=min_g_r(Y,weights,N,q,field_lambda,coupling_lambda,gap_lambda,r,M,nrGapParam,lH,rH,options);
            w(:,r)=wr;
        end
        toc
        matlabpool('close')
    else
        tic
        for r=1:N
            disp(strcat('Minimizing g_r for node r=',int2str(r)))       
            wr=min_g_r(Y,weights,N,q,field_lambda,coupling_lambda,gap_lambda,r,M,nrGapParam,lH,rH,options);
            w(:,r)=wr;
        end
        toc
    end

%Extract the coupling estimates from w.
    JJ=reshape(w(q+1:q+q^2*(N-1),:),q,q,N-1,N);
    Jtemp1=zeros(q,q,N*(N-1)/2);
    Jtemp2=zeros(q,q,N*(N-1)/2);  
    l=1;
    for i=1:(N-1)
         for j=(i+1):N
            Jtemp1(:,:,l)=JJ(:,:,j-1,i); %J_ij as estimated from from g_i.
	    Jtemp2(:,:,l)=JJ(:,:,i,j)'; %J_ij as estimated from from g_j.
            l=l+1;
        end
    end

% Extract G
G=w(q+q^2*(N-1)+1:end,:);


%A note on gauges: 
%The parameter estimates coming from g_r satisfy the gauge
%	lambda_J*sum_s Jtemp_ri(s,k) = 0
%	lambda_J*sum_k Jtemp_ri(s,k) = lambda_h*htemp_r(s)	
%	sum_s htemp_r(s) = 0.
%Only the couplings are used in what follows.
    
    
%Shift the coupling estimates into the Ising gauge.
    J1=zeros(q,q,N*(N-1)/2);
    J2=zeros(q,q,N*(N-1)/2);
    for l=1:(N*(N-1)/2)
        J1(:,:,l)=Jtemp1(:,:,l)-repmat(mean(Jtemp1(:,:,l)),q,1)-repmat(mean(Jtemp1(:,:,l),2),1,q)+mean(mean(Jtemp1(:,:,l)));
	J2(:,:,l)=Jtemp2(:,:,l)-repmat(mean(Jtemp2(:,:,l)),q,1)-repmat(mean(Jtemp2(:,:,l),2),1,q)+mean(mean(Jtemp2(:,:,l)));
    end
%Take J_ij as the average of the estimates from g_i and g_j.
    J=0.5*(J1+J2);

%Calculate frob. norms FN_ij.
    NORMS=zeros(N,N); 
    l=1;
    for i=1:(N-1)
        for j=(i+1):N
            NORMS(i,j)=norm(J(1:end,1:end,l),'fro');
            NORMS(j,i)=NORMS(i,j);
            l=l+1;
        end
    end               
       
    
%Calculate final scores, CN_ij=FN_ij-(FN_i-)(FN_-j)/(FN_--), where '-'
%denotes average.
    norm_means=mean(NORMS)*N/(N-1);
    norm_means_all=mean(mean(NORMS))*N/(N-1);
    CORRNORMS=NORMS-norm_means'*norm_means/norm_means_all;
    output=[];
    for i=1:(N-1)
        for j=(i+1):N
            output=[output;[i,j,CORRNORMS(i,j)]];
        end
    end
    dlmwrite(outputfile,output,'precision',5)
  %  sc=CORRNORMS;
end


function [wr]=min_g_r(Y,weights,N,q,field_lambda,coupling_lambda,gap_lambda,r,M,nrGapParam,lH,rH,options)
%Creates function object for (regularized) g_r and minimizes it using minFunc.
    r=int32(r);
    funObj=@(wr)g_r(wr,Y,weights,N,q,field_lambda,coupling_lambda,gap_lambda,r,M,lH,rH);        
    wr0=zeros(q+q^2*(N-1)+nrGapParam,1);
    wr=minFunc(funObj,wr0,options);    
end

function [fval,grad] = g_r(wr,Y,weights,N,q,lambdah,lambdaJ,lambdaG,r,M,lH,rH)
%Evaluates (regularized) g_r using the mex-file.
	h_r=reshape(wr(1:q),1,q);
	J_r=reshape(wr(q+1:q+q^2*(N-1)),q,q,N-1);
	G=wr((q+q^2*(N-1)+1):end);

	r=int32(r);
	[fval,grad1,grad2,grad3] = g_rC(Y-1,weights,h_r,J_r,[lambdah;lambdaJ;lambdaG],r,G,M,lH,rH);
	grad = [grad1(:);grad2(:);grad3(:)];
end

function [N,B,q,Y] = return_alignment(inputfile)
%Reads alignment from inputfile, removes inserts and converts into numbers.
    align_full = fastaread(inputfile);
    B = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Y = zeros(B,N);

    for i=1:B
        counter = 0;
        for j=1:length(ind)
            if( ind(j) )
                counter = counter + 1;
                Y(i,counter)=letter2number( align_full(i).Sequence(j) );
            end
        end
    end
    q=max(max(Y));
end

function x=letter2number(a)
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'    
            x=2;    
        case 'C'    
            x=3;
        case 'D'
            x=4;
        case 'E'  
            x=5;
        case 'F'
            x=6;
        case 'G'  
            x=7;
        case 'H'
            x=8;
        case 'I'  
            x=9;
        case 'K'
            x=10;
        case 'L'  
            x=11;
        case 'M'
            x=12;
        case 'N'  
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'  
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
end












