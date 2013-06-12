%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%   MutSigCV                                                                          %
%   v1.0                                                                              %
%                                                                                     %
%   (C) 2008-2013 Mike Lawrence & Gaddy Getz                                          %
%   Broad Institute of MIT and Harvard                                                %
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% By downloading the PROGRAM you agree to the following terms of use:
%% 
%% BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
%% FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
%% 
%% This Agreement is made between the Broad Institute, Inc. with a principal 
%% address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the 
%% LICENSEE and is effective at the date the downloading is completed 
%% ("EFFECTIVE DATE").
%% WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, 
%% and BROAD wishes to have this PROGRAM utilized in the public interest, 
%% subject only to the royalty-free, nonexclusive, nontransferable license 
%% rights of the United States Government pursuant to 48 CFR 52.227-14; and
%% WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant 
%% a license on the following terms and conditions.
%% NOW, THEREFORE, in consideration of the promises and covenants made herein, 
%% the parties hereto agree as follows:
%% 
%% 1. DEFINITIONS
%% 1.1	"PROGRAM" shall mean copyright in the object code and source code 
%% known as MutSig and related documentation, if any, as they exist on the 
%% EFFECTIVE DATE and can be downloaded from 
%% http://www.broadinstitute.org/cancer/cga/MutSig on the EFFECTIVE DATE.
%% 
%% 2. LICENSE
%% 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to 
%% LICENSEE, solely for academic non-commercial research purposes, a non-
%% exclusive, non-transferable license to: (a) download, execute and display 
%% the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
%% 
%% LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-
%% free, irrevocable license to any LICENSEE bug fixes or modifications to the 
%% PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE 
%% agrees to provide any such modifications and bug fixes to BROAD promptly 
%% upon their creation.
%% 
%% The LICENSEE may apply the PROGRAM in a pipeline to data owned by users 
%% other than the LICENSEE and provide these users the results of the PROGRAM 
%% provided LICENSEE does so for academic non-commercial purposes only.  For 
%% clarification purposes, academic sponsored research is not a commercial use 
%% under the terms of this Agreement.
%% 
%% 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or 
%% distribute the PROGRAM, in whole or in part, without prior written 
%% permission from BROAD.  LICENSEE shall ensure that all of its users agree 
%% to the terms of this Agreement.  LICENSEE further agrees that it shall not 
%% put the PROGRAM on a network, server, or other similar technology that may 
%% be accessed by anyone other than the LICENSEE and its employees and users 
%% who have agreed to the terms of this agreement.
%% 
%% 2.3  License Limitations. Nothing in this Agreement shall be construed to 
%% confer any rights upon LICENSEE by implication, estoppel, or otherwise to 
%% any computer software, trademark, intellectual property, or patent rights 
%% of BROAD, or of any other entity, except as expressly granted herein. 
%% LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for 
%% any commercial purpose, including without limitation, as the basis of a 
%% commercial software or hardware product or to provide services. LICENSEE 
%% further agrees that the PROGRAM shall not be copied or otherwise adapted in 
%% order to circumvent the need for obtaining a license for use of the 
%% PROGRAM.  
%% 
%% 3. OWNERSHIP OF INTELLECTUAL PROPERTY 
%% LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. 
%% The PROGRAM is marked with the following BROAD copyright notice and notice 
%% of attribution to contributors. LICENSEE shall retain such notice on all 
%% copies.  LICENSEE agrees to include appropriate attribution if any results 
%% obtained from use of the PROGRAM are included in any publication.
%% 
%% Copyright 2012 Broad Institute, Inc.
%% Notice of attribution:  The MutSig program was made available through the 
%% generosity of the Cancer Genome Analysis group at the Broad Institute, Inc. 
%% 
%% LICENSEE shall not use any trademark or trade name of BROAD, or any 
%% variation, adaptation, or abbreviation, of such marks or trade names, or 
%% any names of officers, faculty, students, employees, or agents of BROAD 
%% except as states above for attribution purposes.
%% 
%% 4. INDEMNIFICATION
%% LICENSEE shall indemnify, defend, and hold harmless BROAD, and their 
%% respective officers, faculty, students, employees, associated investigators 
%% and agents, and their respective successors, heirs and assigns, 
%% ("Indemnitees"), against any liability, damage, loss, or expense (including 
%% reasonable attorneys fees and expenses) incurred by or imposed upon any of 
%% the Indemnitees in connection with any claims, suits, actions, demands or 
%% judgments arising out of any theory of liability (including, without 
%% limitation, actions in the form of tort, warranty, or strict liability and 
%% regardless of whether such action has any factual basis) pursuant to any 
%% right or license granted under this Agreement.
%% 
%% 5. NO REPRESENTATIONS OR WARRANTIES
%% THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR 
%% WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR 
%% IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, 
%% FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT 
%% OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES 
%% OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER 
%% LITERATURE MAY BE ISSUED FROM TIME TO TIME.
%% IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, 
%% AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR 
%% CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC 
%% DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD 
%% SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF 
%% THE POSSIBILITY OF THE FOREGOING.
%% 
%% 6. ASSIGNMENT
%% This Agreement is personal to LICENSEE and any rights or obligations 
%% assigned by LICENSEE without the prior written consent of BROAD shall be 
%% null and void.
%% 
%% 7. MISCELLANEOUS
%% 7.1 Export Control. LICENSEE gives assurance that it will comply with all 
%% United States export control laws and regulations controlling the export of 
%% the PROGRAM, including, without limitation, all Export Administration 
%% Regulations of the United States Department of Commerce. Among other 
%% things, these laws and regulations prohibit, or require a license for, the 
%% export of certain types of software to specified countries.
%% 7.2 Termination. LICENSEE shall have the right to terminate this Agreement 
%% for any reason upon prior written notice to BROAD. If LICENSEE breaches any 
%% provision hereunder, and fails to cure such breach within thirty (30) days, 
%% BROAD may terminate this Agreement immediately. Upon termination, LICENSEE 
%% shall provide BROAD with written assurance that the original and all copies 
%% of the PROGRAM have been destroyed, except that, upon prior written 
%% authorization from BROAD, LICENSEE may retain a copy for archive purposes.
%% 7.3 Survival. The following provisions shall survive the expiration or 
%% termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 
%% 7.3, and 7.4.
%% 7.4 Notice. Any notices under this Agreement shall be in writing, shall 
%% specifically refer to this Agreement, and shall be sent by hand, recognized 
%% national overnight courier, confirmed facsimile transmission, confirmed 
%% electronic mail, or registered or certified mail, postage prepaid, return 
%% receipt requested.  All notices under this Agreement shall be deemed 
%% effective upon receipt. 
%% 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, 
%% supplemented, or otherwise modified only by means of a written instrument 
%% signed by all parties. Any waiver of any rights or failure to act in a 
%% specific instance shall relate only to such instance and shall not be 
%% construed as an agreement to waive any rights or fail to act in any other 
%% instance, whether or not similar. This Agreement constitutes the entire 
%% agreement among the parties with respect to its subject matter and 
%% supersedes prior agreements or understandings between the parties relating 
%% to its subject matter. 
%% 7.6 Binding Effect; Headings. This Agreement shall be binding upon and 
%% inure to the benefit of the parties and their respective permitted 
%% successors and assigns. All headings are for convenience only and shall not 
%% affect the meaning of any provision of this Agreement.
%% 7.7 Governing Law. This Agreement shall be construed, governed, interpreted 
%% and applied in accordance with the internal laws of the Commonwealth of 
%% Massachusetts, U.S.A., without regard to conflict of laws principles.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MutSigCV(mutation_file,coverage_file,covariate_file,output_file)

if nargin~=4, error('usage: MutSigCV(mutation_file,coverage_file,covariate_file,output_file)'); end

% first, ensure output_file directory exists
[outpath name ext versn] = fileparts(output_file);
if ~isempty(outpath) && ~exist(outpath,'dir'), mkdir(outpath); end

fprintf('Loading data...\n');

M = load_struct(mutation_file);
M = make_numeric(M,'categ');
M = make_boolean(M,{'is_coding','is_silent'});
G = load_struct_specify_string_cols(covariate_file,1);
C = load_struct_specify_string_cols(coverage_file,1:2);

fprintf('Building n and N tables...\n');

f = fieldnames(C); patient_names = f(4:end);
f = fieldnames(G); cvnames = f(2:end);

ng = slength(G);
np = length(patient_names);
ncat = max(M.categ);
nv = length(cvnames);

% make sure C is sorted by the same gene order as in G
C.gene_idx = listmap(C.gene,G.gene);
C = sort_struct(C,'gene_idx');

M.gene_idx = listmap(M.Hugo_Symbol,G.gene);
M.patient = regexprep(M.Tumor_Sample_Barcode,'-Tumor$','');
M.patient = regexprep(M.patient,'-','_');
M.patient_idx = listmap(M.patient,patient_names);

midx = find(M.is_coding & M.is_silent);
n_silent = hist3d(M.gene_idx(midx),M.categ(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
midx = find(M.is_coding & ~M.is_silent);
n_nonsilent = hist3d(M.gene_idx(midx),M.categ(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
midx = find(~M.is_coding);
n_flank = hist3d(M.gene_idx(midx),M.categ(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);

N_silent = nan(ng,ncat,np);
N_nonsilent = nan(ng,ncat,np);
N_flank = nan(ng,ncat,np);

for ci=1:ncat
  silent_idx = find(strcmp(C.zone,'silent') & C.categ==ci);
  nonsilent_idx = find(strcmp(C.zone,'nonsilent') & C.categ==ci);
  flank_idx = find(strcmp(C.zone,'flank') & C.categ==ci);
  for pi=1:np
    N_silent(:,ci,pi) = C.(patient_names{pi})(silent_idx);
    N_nonsilent(:,ci,pi) = C.(patient_names{pi})(nonsilent_idx);
    N_flank(:,ci,pi) = C.(patient_names{pi})(flank_idx);
  end
end

% add total columns
n_silent(:,end+1,:) = sum(n_silent,2);
n_nonsilent(:,end+1,:) = sum(n_nonsilent,2);
n_flank(:,end+1,:) = sum(n_flank,2);
N_silent(:,end+1,:) = N_silent(:,end,:);          % copy total coverage from null+indel coverage
N_nonsilent(:,end+1,:) = N_nonsilent(:,end,:);
N_flank(:,end+1,:) = N_flank(:,end,:);

% total across patients, save in G
G.N_nonsilent = sum(N_nonsilent(:,end,:),3);
G.N_silent = sum(N_silent(:,end,:),3);
G.N_flank = sum(N_flank(:,end,:),3);
G.n_nonsilent = sum(n_nonsilent(:,end,:),3);
G.n_silent = sum(n_silent(:,end,:),3);
G.n_flank = sum(n_flank(:,end,:),3);

fprintf('Processing covariates...\n');

V = nan(ng,nv);
for vi=1:nv, V(:,vi) = G.(cvnames{vi}); end

% convert covariate raw values to Z-scores
Z = nan(ng,nv);
for vi=1:nv
  missing = isnan(V(:,vi)) | isinf(V(:,vi));
  mn = mean(V(~missing,vi));
  sd = std(V(~missing,vi));
  Z(~missing,vi) = (V(~missing,vi)-mn)./sd;
end

fprintf('Finding bagels...  ');

max_neighbors = 50;
qual_min = 0.05;

G.nnei = nan(ng,1); G.x = nan(ng,1); G.X = nan(ng,1);

for gi=1:ng, if ~mod(gi,1000), fprintf('%d/%d ',gi,ng); end

  % calculate distances from this gene
  df2 = bsxfun(@minus,Z,Z(gi,:)).^2;
  dist2 = nansum(df2,2)./sum(~isnan(df2),2);
  [tmp,ord] = sort(dist2); ord = [gi;ord(ord~=gi)];

  % expand bagel outward until quality falls below qual_min
  nfit=0; Nfit=0;
  for ni=0:max_neighbors, gidx = ord(ni+1);

    ngene = G.n_silent(gidx) + G.n_flank(gidx);
    Ngene = G.N_silent(gidx) + G.N_flank(gidx);
    if ni==0, ngene0=ngene; Ngene0=Ngene; end
    nfit=nfit+ngene; Nfit=Nfit+Ngene;

    % compare the gene being added to the central gene
    hc = hyge2cdf(ngene,Ngene,ngene0,Ngene0);
    qual_left = min(hc, 1-hc);
    qual = 2*qual_left;

    % stopping criterion: stop if this gene would drop quality below qual_min
    if ni>0 && qual<qual_min, break; end

    % update gene's statistics
    G.nnei(gi) = ni; G.x(gi) = nfit; G.X(gi) = Nfit;

  end % next neighborhood size
end, fprintf('\n'); % next gene

fprintf('Expanding to (x,X)_gcp...\n');

n_gcp = n_nonsilent + n_silent + n_flank;
N_gcp = N_nonsilent + N_silent + N_flank;

n_cp = sum(n_gcp,1);
N_cp = sum(N_gcp,1);

n_c = sum(n_cp,3);
N_c = sum(N_cp,3);
mu_c = n_c./N_c;

n_tot = n_c(end);
N_tot = N_c(end);
mu_tot = n_tot/N_tot;
f_c = mu_c/mu_tot;
f_Nc = N_c/N_tot;

n_p = n_cp(:,end,:);
N_p = N_cp(:,end,:);
mu_p = n_p./N_p;
f_p = mu_p/mu_tot;
f_Np = N_p/mean(N_p);

x_gcp = repmat(G.x,[1 ncat+1 np]); X_gcp = repmat(G.X,[1 ncat+1 np]);       % last column = total
x_gcp = bsxfun(@times,x_gcp,f_c.*f_Nc); X_gcp = bsxfun(@times,X_gcp,f_Nc);
x_gcp = bsxfun(@times,x_gcp,f_p.*f_Np); X_gcp = bsxfun(@times,X_gcp,f_Np);

fprintf('Calculating p-value using 2D Projection method...  ');

null_score_boost = 3;
min_effect_size = 1.25;
convolution_numbins = 1000;

G.p = nan(ng,1);

for g=1:ng, if ~mod(g,1000), fprintf('%d/%d ',g,ng); end

  % STEP 1
  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.

  N = reshape(N_nonsilent(g,1:ncat,:),ncat,np)';
  n = reshape(n_nonsilent(g,1:ncat,:),ncat,np)';
  x = reshape(x_gcp(g,1:ncat,:),ncat,np)';
  X = reshape(X_gcp(g,1:ncat,:),ncat,np)';
  P0 = hyge2pdf(0,N,x,X);
  P1 = hyge2pdf(1,N,x,X);

  % determine each patient's priority order of categories (according to P1)
  % left column of "priority" = least extreme category of mutation
  % right column of "priority" = most extreme category of mutation
  [tmp priority] = sort(P1,2,'descend');
  % sort the P arrays to put the columns in least->most priority order
  shft = (priority - repmat(1:ncat,np,1));
  map = reshape(1:(np*ncat),np,ncat);
  newmap = map + shft*np;
  P0 = P0(newmap);
  P1 = P1(newmap);
  P2 = 1-(P0+P1);  % note, P2 means "P(2+)"
  P2(P2<0) = 0;

  % STEP 2
  % for each sample, compute probability that it would have been of each (2-dimensional) degree.
  % degree=(d1,d2), where d=0 (no mut) ..... ncat (most extreme mut)
  % d1 is the MOST extreme mutation (or no mutation)
  % d2 is the SECOND MOST extreme mutation (or no mutation)
  % d1 can be 0-ncat; d2 can be 0-d1

  Pdeg = zeros(np,ncat+1,ncat+1);
  for d1=0:ncat, for d2=0:d1
    % has to have 0 in any/all categories > d1
    p = prod(P0(:,d1+1:end),2);
    if (d1>0)  % and (if d1>0)
      if (d1==d2)
        % if d1==d2, has to have 2+ in category d1
        p = p .* P2(:,d1);
      else
        % else:   has to have exactly 1 in category d1
        %         has to be clear in any/all categories (d2+1) to (d1-1)
        %         and (if d2>0) have (1 or 2+) in category d2
        p = p .* P1(:,d1);
        p = p .* prod(P0(:,d2+1:d1-1),2);
        if (d2>0)
          p = p .* (P1(:,d2)+P2(:,d2));
        end
      end
    end
    Pdeg(:,d1+1,d2+1) = p;
  end,end

  %% STEP 2a: calculate score for a sample being of each possible degree
  %% (uses new style, where score = -log10 probability of having a mutation in that category
  %% (zero score for no mutations)
  Sdeg = zeros(np,ncat+1,ncat+1);
  for d1=1:ncat, for d2=0:d1
    if d1==d2
      p = P2(:,d1);
    else
      if d2>0
        p = P1(:,d1).*P1(:,d2);
      else
        p = P1(:,d1);
      end
    end
    Sdeg(:,d1+1,d2+1) = -log10(p);
  end,end

  % null score boost
  priority2 = [zeros(np,1) priority];
  Sdeg(priority2==ncat) = Sdeg(priority2==ncat) + null_score_boost;

 % STEP 3
  % determine actual (two-dimensional) degree and score for each sample
  % sum scores to get score_obs for gene

  degree = zeros(np,2);
  score_obs = 0;
  for p = 1:np
    i = 1;
    for d = ncat:-1:1
      c = priority(p,d);
      if i==1
        if n(p,c)>=2
          degree(p,:) = [d d];
          i = 3;
        elseif n(p,c)==1
          degree(p,i) = d;
          i=i+1;
        end
      elseif i==2
        if n(p,c)>=1
          degree(p,i) = d;
          i=i+1;
        end
      else % i>2: done
        break
      end
    end
    score_sample = Sdeg(p,degree(p,1)+1,degree(p,2)+1);
    score_obs = score_obs + score_sample;
  end

  % minimum effect size
  score_obs = score_obs / min_effect_size;

  % for zero score, don't bother doing convolutions
  if score_obs<=0, G.p(g)=1; continue; end

  % STEP 4
  % compute P value for gene by convolutions

  numbins = convolution_numbins;
  binsize = score_obs / numbins;
  H = zeros(numbins,1);
  H(1) = 1;  % initial condition: all probability is in first bin

  % sequential convolution
  offset = min(numbins, round(Sdeg/binsize));
  ncols = (ncat+1)*(ncat+2)/2;
  newH = zeros(numbins,ncols);
  for p=1:np
    newH(:) = 0;
    col=1;
    for d1=0:ncat, for d2=0:d1
      o = offset(p,d1+1,d2+1);
      newH(o+1:end,col) = Pdeg(p,d1+1,d2+1) .* H(1:end-o);
      col=col+1;
    end,end
    H = sum(newH,2);
  end

  % save p-value
  G.p(g) = max(0,1-sum(H));

end, fprintf('\n');   % next gene

% FDR
G.q = calc_fdr_value(G.p);

G = sort_struct(G,'p');
save_struct(G,output_file);

fprintf('Done.\n');

end

% end of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% math subfunctions                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=hyge2pdf(k,n,k1,n1)
  p = exp(gammaln(n1+2) - gammaln(k1+1) - gammaln(n1-k1+1) + ...
          gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
          gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2));
  p = max(0,min(1,p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=hyge2cdf(k,n,k1,n1)
  p=0;
  for ki=0:k, p=p+hyge2pdf(ki,n,k1,n1); end
  p = max(0,min(1,p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fdr=calc_fdr_value(p)
  if isempty(p)
    fdr = p;
    return
  end
  if size(p,1)==1
    p=p';
    trans=1;
  else
    trans=0;
  end
  [sp,ord]=sort(p);
  fdr=sp*size(p,1)./repmat((1:size(p,1))',1,size(p,2));
  fdr(fdr>1)=1;
  fdr=[fdr; ones(1,size(fdr,2))];
  for i=size(p,1):-1:1
    fdr(i,:)=min(fdr(i:(i+1),:),[],1);
  end
  fdr=fdr(1:(end-1),:);
  ordmat=ord+repmat(0:size(p,1):size(p,1)*(size(fdr,2)-1),size(ord,1),1);
  fdr(ordmat(:))=fdr(:);
  if trans
    fdr=fdr';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = hist3d(a,b,c,firsta,lasta,firstb,lastb,firstc,lastc)
  if nargin~=9, error('requires 9 input arguments'); end
  h = zeros(lasta-firsta+1,lastb-firstb+1,lastc-firstc+1);
  for bi=firstb:lastb, for ci=firstc:lastc
      h(:,bi,ci) = histc(a(b==bi & c==ci),firsta:lasta);
  end,end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% technical subfunctions: file input/output, etc.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ok = demand_file(fname)
  if ~iscell(fname), fname = {fname}; end
  ok = true(size(fname));
  for i=1:numel(fname), if ~mod(i,100), fprintf('%d/%d ',i,numel(fname)); end
    if ~exist(fname{i},'file'), ok(i)=false; end
  end, if numel(fname)>=100, fprintf('\n'); end
  nnf = sum(~ok);
  nf = sum(ok);
  if nargout==0
    if nf>=10, fprintf('%d files found\n',nf); end
  end
  if nnf>0
    blank = 0;
    for i=1:numel(fname)
      if ~ok(i)
        if strcmp(fname{i},'')
          blank=blank+1;
        else
          if nnf<20
            fprintf('\tNot found: %s\n',fname{i});
          end
        end
      end
    end
    if blank>0, fprintf('\tNot found: <blank> (%d)\n',blank); end
    if nargout==0
      if nnf==1, error('1 file not found'); else error('%d files not found',nnf); end
    end
  end
  if nargout==0, clear ok; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=dlmsep(s,d)
  if nargin==1
    d=9; % tab
  end

  pos=find(ismember(s,d));
  if ~isempty(pos)
    pos=[ 0 pos length(s)+1];
    for i=1:(length(pos)-1)
      res{i}=s((pos(i)+1):(pos(i+1)-1));
    end
  else
    res{1}=s;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = ensure_fopen(varargin)
  result = fopen(varargin{:});
  if (result==-1), error('ERROR WRITING TO %s',varargin{1}); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d=find_dlm(fname,dlm)
  if ~exist('dlm','var') || isempty(dlm)
    dlm=[ char(9) ',|'];
  end
  f=fopen(fname,'r');
  l=fgetl(f);
  for i=1:length(dlm)
    h(i)=length(find(l==dlm(i)));
  end
  [hm,hi]=max(h);
  d=dlm(hi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = flip(A)
  if sum(size(A)>1)>1, error('flip only works for vectors'); end
  if size(A,1)==1, A=fliplr(A); else A=flipud(A); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function numcols = get_colcount(fname,num_header_lines)
  if ~exist('num_header_lines','var'), num_header_lines = 1; end
  if ~exist(fname,'file'), error('%s not found',fname); end
  
  f = fopen(fname);
  for i=1:num_header_lines+1; l = fgetl(f); end
  numcols = sum(l==char(9))+1;
  fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=as_column(x)
  if size(x,2)>1
    x=x';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = impose_default_value(S,field,value,acceptable_values)
  if ~isfield(S,field) || isempty(getfield(S,field))
    if ischar(value) & strcmp(value,'*required*')
      error('%s is a required field of P',field);
    else
      try
        S=setfield(S,field,value);
      catch me
        fprintf('Error setting field "%s"\n',field);
        disp(me);disp(me.message);
      end
    end
  end
  
  if exist('acceptable_values','var')
    av = acceptable_values;
    v = getfield(S,field);
    if ischar(v) && isnumeric(av)
      v = str2double(v);
      S = setfield(S,field,v);
    end
    if ~ischar(v) && ischar(av)
      error('%s is assigned value of different type than acceptable_values',field);
    end
    try
      ism = ismember(v,av);
    catch me
      error('%s is assigned value of different type than acceptable_values',field);
    end
    if ~ism
      fprintf('Acceptable values for %s:\n',field); disp(av);
      fprintf('Attempted to set to:\n'); disp(v)
      error('Invalid setting of %s',field);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = listmap(a,b)
  if ischar(a), a={a}; end
  if ischar(b), b={b}; end
  m = nan(length(a),1);
  [a ai aj] = unique(a);
  [c ia ib] = intersect(a,b);
  for i=1:length(c), if ~mod(i,1e5), fprintf('%d/%d ',i,length(c)); end
    m(aj==ia(i)) = ib(i);
  end, if length(c)>=1e5, fprintf('\n'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = load_lines(fname)
  X = load_textfile(fname);
  L = text_to_lines(X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = load_struct(varargin)
% load_struct(filename, format, separator_character, header_lines, <other_parameters>, P)
  
% loads a tab-delimited table from the specified file
% into a structure with field names based on the column headings
%
% piggybacks on read_table()
% all parameters are passed to read_table(),
% except last parameter if it is a struct (i.e. P struct)
% P struct can have following members:
%   P.lowercase_fieldnames: if true or 1, fieldnames are converted to lowercase
%
% WAS: load_struct(filename, format, header_lines, lowercase_fieldnames)
%   format is [] by default.
%   header_lines is 1 by default.
%   if header_lines is set to 0, field names are col1, col2, etc.
%   if lowercase_fieldnames is 1, fieldnames are converted to lowercase
%
% Mike Lawrence 2008-04-29
% fixed 2010-12-10 to give all its parameters to read_table
%   (except last parameter if a P struct)

% check parameters

  args = varargin;
  if length(args)<1
    error('filename required');
  end
  filename = args{1};
  if ~ischar(filename)
    error('first parameter should be filename (character string)');
  end
  demand_file(filename);
  
  % see if last parameter is a P struct
  if isstruct(args{end})
    P = args{end};
    args = args(1:end-1);
  else
    P = [];
  end
  
  P = impose_default_value(P,'lowercase_fieldnames',false);
  P = impose_default_value(P,'ignore_poundsign_lines',true);
  
  if length(args)>=2 
    % ok
  end
  if length(args)>=3
    if isempty(args{3})
      error('third parameter should be separator character, e.g. char(9)');
    end
    if isnumeric(args{3})
      error('header_lines is now the fourth parameter... please update call to load_struct');
    end
  end
  if length(args)>=4
    if islogical(args{4})
      error('lowercase_fieldnames has been moved to P struct... please update call to load_struct');
    end
    if ~isnumeric(args{4})
      error('fourth parameter should be number of header lines');
    end
  end
  
  % HANDLE COMMENT LINES:
  default_header_lines = 1;
  %% see if table has comment lines at the beginning (start with #): if so, increment header_lines to skip them
  n_comment_lines = 0;
  
  if P.ignore_poundsign_lines
    f = fopen(args{1});
    while(true)
      x = fgetl(f);	
      if (isempty(x))||(x(1)=='#') % skip empty lines or lines starting with #
        n_comment_lines = n_comment_lines + 1;
        continue
      elseif strncmp(x,'Oncotator v',11)
        fprintf('Un-poundsigned Oncotator header detected and skipped.\n');
        n_comment_lines = n_comment_lines + 1;
        continue
      else
        break
      end
    end
    fclose(f);
  end
  
  default_header_lines = 1 + n_comment_lines;
  
  %% default args
  if length(args)==1, args = [args {''}]; end         % format string
  if length(args)==2, args = [args {char(9)}]; end       % separator character
  if length(args)==3, args = [args {default_header_lines}]; end          % number of header lines
  
  % default whitespace
  has_already = false;
  for i=4:length(args)
    if ischar(args{i}) & strcmpi(args{i},'whitespace'), has_already=true; break; end
  end
  if ~has_already, args = [args {'whitespace'} {'\b\r'}]; end
  
  % default bufSize
  has_already = false;
  for i=4:length(args)
    if ischar(args{i}) & strcmpi(args{i},'bufSize'), has_already=true; break; end
  end
  if ~has_already, args = [args {'bufSize'} {50000}]; end
  
  % LOAD TABLE
  try
    table = read_table(args{:});
    nf = length(table.dat);
  catch me
    q = load_lines(args{1});
    if isempty(q)
      fprintf('\n%s is a blank file\n',args{1});
      table = [];
      table.dlm  = args{3};
      table.headers = {{}};
      table.dat = {{}};
      nf = 0;
    else
      disp(me);
      disp(me.message);
      error('Error loading struct file');
    end
  end
  
  if isempty(table.headers)
    table.headers{1} = cell(nf,1);
    for f=1:nf
      table.headers{1}{f} = sprintf('col%d', f);
    end
  end
  
  % process header line
  
  fields = table.headers{end};
  if length(fields)~=nf
    fprintf('Header line has %d column headers instead of the expected %d:\n',length(fields),nf);
    fields{:}
    error('Unable to parse table header line.');
  end
  
  % remove illegal characters from column headings
  % and convert to list of unique field names
  
  if P.lowercase_fieldnames, fields = lower(fields); end
  fields = regexprep(fields, '\W','');   % remove any characters except A-Z, a-z, 0-9, underscore
  fields_orig = fields;
  fields = genvarname(fields_orig);
  
  % preserve "end", because it's only going to be a field name, not a variable name
  for f=1:nf
    if strcmp(fields_orig{f}, 'end')
      fields{f} = 'end';
      break
    end
    if strcmp(fields_orig{f}, 'End')
      if P.lowercase_fieldnames
        fields{f} = 'end';
      else
        fields{f} = 'End';
      end
      break
    end
  end
  
  S = struct();
  for f=1:nf
    S = setfield(S, fields{f}, table.dat{f});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = load_struct_specify_string_cols(fname,string_cols,num_header_lines,lowercase_fieldnames)

  if ~exist('fname','var'), error('Must specify fname'); end
  if ~exist('string_cols','var'), string_cols = []; end
  if ~exist('num_header_lines','var'), num_header_lines = 1; end
  if ~exist('lowercase_fieldnames','var'), lowercase_fieldnames = false; end
  
  if ~exist(fname,'file'), error('%s not found',fname); end
  
  numcols = get_colcount(fname,num_header_lines);
  
  is_string = false(1,numcols);
  is_string(string_cols) = true;
  format = [];
  for i=1:numcols
    if is_string(i), format = [format '%s'];
    else format = [format '%f']; end
  end
  
  P=[]; P.lowercase_fieldnames = lowercase_fieldnames;
  X = load_struct(fname,format,char(9),num_header_lines,P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t = load_textfile(filename)
  in = fopen(filename);
  t = fread(in,'uint8=>char')';
  fclose(in);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = make_boolean(S,varargin)
  fields = {};
  for i=1:length(varargin)
    fields = [fields varargin{i}];
  end
  for i=1:length(fields)
    if isfield(S,fields{i})    
      x = getfield(S,fields{i});
      if ~islogical(x) && ~isnumeric(x)
        x = str2double(x);
      end
      if ~islogical(x)
        x = (x~=0);
      end
      S = setfield(S,fields{i},x);
    else
      fprintf('No such field: %s\n', fields{i});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = make_numeric(S,varargin)
  fields = {};
  for i=1:length(varargin)
    fields = [fields varargin{i}];
  end
  for i=1:length(fields)
    if isfield(S,fields{i})    
      x = getfield(S,fields{i});
      if isnumeric(x) || islogical(x)
      else
        x = str2double(x);
        S = setfield(S,fields{i},x);
      end
    else
      fprintf('No such field: %s\n', fields{i});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = nansub(X,idx,filler)
  if iscellstr(X) && size(X,1)==1 && size(X,2)>1
    fprintf('WARNING: current implementation has problem with row-vector cellstr arrays\n');
  end
  if islogical(X)
    type = 0;
  elseif isnumeric(X)
    type = 1;
  elseif iscell(X)
    type = 2;
  else
    error('Unsuuported array type');
  end
  
  if ~exist('filler','var')
    if type==0
      filler = false;
    elseif type==1
      filler = nan;
    elseif type==2
      filler = {''};
    else
      error('Inconsistent behavior with "type"');
    end
  end
  
  if type==0
    if ~islogical(filler)
      error('Inappropriate filler for logical array');
    end
  elseif type==1
    if ~isnumeric(filler)
      error('Inappropriate filler for numeric array');
    end
  elseif type==2
    if ischar(filler)
      filler = {filler};
    end
    if ~iscell(filler)
      error('Inappropriate filler for cell array');
    end
  else
    error('Inconsistent behavior with "type"');
  end
  
  sz = size(X); sz(1) = length(idx);
  Y = repmat(filler,sz);
  idx2 = find(~isnan(idx) & idx>=1 & idx<=length(X));
  Y(idx2,:,:,:,:) = X(idx(idx2),:,:,:,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function st=newline
  if ispc
    st='\r\n';
  else
    st='\n';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = parse(S,r,f,numeric)
  if ~iscell(S), S = tolines(S); end
  tokens = regexp(S,r,'tokens');
  if ~iscell(f), f = {f}; end
  x = tokens2struct(tokens,f);
  if exist('numeric','var'), x = make_numeric(x,f(numeric)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fl,fid]=read_dlm_file(fname,dlm,nlines)
% READ_DLM_FILE  Read all or some of a delimited file into cell arrays-of-arrays.
% [FL,FID] = READ_DLM_FILE(FNAME,DLM,NLINES)
%    Reads a file (FNAME) and partitions to line using the delimeter (DLM).
%    The default DLM value is \t (9). The optional NLINES parameter
%    limits the number of lines read. 
%    
%    FL is a cell array of cell arrays of strings.  FL{n} is a cell array 
%    strings containing the nth line of the file.  
%    FID is the fid for the opened file, FNAME.
%
%
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%

  if nargin==1
    dlm=9;
  end
  
  if ischar(fname)
    fid=fopen(fname,'r');
  else
    fid=fname;
  end
  
  ln=1;
  if ~exist('nlines','var')
    nlines=Inf;
    do_close=1;
  else
    do_close=0;
  end
  
  had_output=0;
  fl={};
  while(ln<=nlines)
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fl{ln}=dlmsep(tline,dlm);         
    ln=ln+1;
    if mod(ln,1000)==0
      verbose(['...' num2str(ln)],30);
      had_output=1;
    end
  end
  if do_close
    fclose(fid);
    fid=-1;
  end
  ln=ln-1;
  if had_output
    verbose([newline],30);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tab=read_table(fname,format,dlm,headerlines,varargin)

  fpos=0;
  if headerlines~=0
    if ischar(fname)
      f=fopen(fname,'r');
    else
      f=fname;
      fpos=ftell(f);
    end
    if length(dlm)~=1
      tab.dlm=find_dlm(fname,dlm);
    else
      tab.dlm=dlm;
    end
    if headerlines>0
      tab.headers=read_dlm_file(f,dlm,headerlines);
    elseif headerlines==-1 % R convention
      headerlines=1;
      tab.headers=read_dlm_file(f,dlm,headerlines);
      tab.headers{1}=['EMPTY' tab.headers{1,:}];
    end
  else
    if ischar(fname)
      f=fopen(fname,'r');
    else
      f=fname;
      fpos=ftell(f);
    end
    tab.headers={};
    tab.dlm=dlm;
  end
  
  if isempty(format)
    if isempty(tab.headers)
      error('must have either format or headerlines');
    else
      format=[repmat('%s',1,length(tab.headers{end})) '\n'];   % (to allow for multiple header lines)
    end
  elseif iscell(format)
    if isempty(tab.headers)
      error('must have either format or headerlines');
    else
      if length(format)==1
        format=[repmat(format{1},1,length(tab.headers{end})) '\n'];
      else
        format=[format{1} repmat(format{3},1,length(tab.headers{end})-format{2}) '\n'];
      end
    end
  end

  if strcmp(format((end-1):end),'\n')
    format=format(1:(end-2));
  end
  
  verbose(['Reading file using format:' format],10);
  fseek(f,fpos,'bof');
  tab.dat=textscan(f,format,'headerLines',headerlines,'delimiter',tab.dlm,varargin{:});
  fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s order]=reorder_struct(s,order)
  
  if nargin~=2, error('reorder_struct(s,order)'); end
  
  if islogical(order), order = find(order); end
  if ischar(order)
    if strcmpi(order,'end')
      order = slength(s);
    else
      error('invalid index parameter');
    end
  end
  
  order = as_column(order);
  nanflag = any(isnan(order));
  fields = fieldnames(s);
  nf = length(fields);

  for i=1:nf
    f = getfield(s,fields{i});
    if nanflag
      f = nansub(f,order);
    else
      f = f(order,:,:,:,:,:,:,:,:,:);
    end
    s = setfield(s,fields{i},f);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = save_struct(S,filename,noheader,deprecated_formatflag)
%
% save_struct(S, filename)
%
% writes a tab-delimited table from the given struct.
% writes a header line based on the field names.
%
% fields are written in the order in which they appear in the struct.

  if nargin==0, error('requires an input argument'); end
  if nargin>4, error('too many input arguments'); end
  if nargin==4, fprintf('USE OF DEPRECATED FORMATFLAG\n'); end
  
  if nargout==0 && nargin<2, error('requires at least 2 arguments'); end
  if nargout==1 && (nargin<2 || isempty(filename))
    return_string_only = true;
  else
    return_string_only = false;
  end
  if nargout>1, error('only a single output possible'); end
  
  if exist('noheader','var') && ~isempty(noheader) && (~islogical(noheader) || ~(noheader==false))
    if (islogical(noheader) && noheader==true) || strncmpi(noheader,'no_header',9) || strncmpi(noheader,'noheader',8)
      noheader = true;
    else
      error('third parameter should be "no_headers", true, false, empty, or nothing.');
    end
  else
    noheader = false;
  end
  
  if ~return_string_only
    out = ensure_fopen(filename,'wt');
  end
  
  if ~isstruct(S)
    if isempty(S)
      S = struct;
    else
      error('S should be a struct');
    end
  end
  
  fld = fieldnames(S);
  nf = length(fld);
  slen = slength(S);
  
  % see if struct is empty
  if slen==0
    if return_string_only
      F = '';
    else
      for i=1:nf
        fprintf(out,fld{i});
        if i==nf
          fprintf(out,'\n');
        else
          fprintf(out,'\t');
        end
      end
      fclose(out);
    end
    return
  end

  % see if struct is to long too handle all at once
  chunksize = round(1e7/nf);
  if slen>chunksize
    if return_string_only
      error('struct is too large to convert all at once in memory');
    end
    
    for st=1:chunksize:slen
      fprintf('HUGE STRUCT: SAVING CHUNK %d/%d\n', ceil(st/chunksize), ceil(slen/chunksize));
      en = min(slen,st+chunksize-1);
      Si = reorder_struct(S,st:en);
      if exist('deprecated_formatflag','var')
        F = save_struct(Si,[],noheader,deprecated_formatflag);
      else
        F = save_struct(Si,[],noheader);
      end
      fwrite(out,F);
      clear F;
      noheader = 'no_header'; % subsequent chunks omit header
    end
    
  else   % struct is not too big to save all at once
    F = cell(1,nf*2);
    nr = -1;
    tt0 = tic;
      
      for f=1:nf
        tt1 = toc(tt0);
        if tt1>10
          if ~exist('flag01','var')
            fprintf('  [save_struct] ');
            flag01 = true;
          end
          fprintf('%d/%d ',f,nf);
        end
        C = getfield(S,fld{f});    % get column
                                   % check for legal type
        if isempty(C)
          C = {};
        elseif isnumeric(C) || islogical(C)
          if ndims(C)>2
            fprintf('Field %s is multidimensional: skipping\n', fld{f});
            C = {};
          elseif size(C,2)>1
            fprintf('Field %s is not a column vector: skipping\n', fld{f});
            C = {};
          else
            if exist('deprecated_formatflag','var') && deprecated_formatflag     % convert to cell array of strings
              C = cellstr(num2str(C));
            else                       % note: "-" is important to avoid extra spaces in output
              if any(mod(C,1))
                C = cellstr(num2str(C,'%-d'));     % not all integers
              else
                C = cellstr(num2str(C,'%-.0f'));    % all integers
              end
            end
          end
        else  % column is a cell
          idx = find(cellfun(@iscell,C));
          if ~isempty(idx)
            fprintf('WARNING: Field %s contains %d entries that are cells:\n', fld{f}, length(idx));
            fprintf('         replacing these with "?????{cell}?????"\n');
            C(idx) = repmat({'?????{cell}?????'},length(idx),1);
          end
          idx = find(~cellfun(@isempty,C) & ~cellfun(@ischar,C));
          if ~isempty(idx)
            fprintf('WARNING: Field %s contains %d entries that are not chars:\n', fld{f}, length(idx));
            fprintf('         replacing these with "?????{unknown_format}?????"\n');
            C(idx) = repmat({'?????{unknown_format}?????'},length(idx),1);
          end
        end
        if isempty(C)
          if nr==-1, error('Problematic leftmost column'); end
          C = repmat({'...'},nr,1);
        end
        % check for consistent column length
        if nr==-1, nr=length(C); end
        if nr~=length(C), error('Field %s is a different length', fld{f}); end
        % add column title
        if ~noheader, C = [fld{f}; C]; end
        % add column to file
        F{f*2-1} = C;
        % add tab or newline
        F{f*2} = repmat({char(9+(f==nf))},length(F{f*2-1}),1);
      end
      
      if tt1>10, fprintf(' [collapse]'); end
      F = strcat(F{:});              % collapse columns to lines
      F = [F{:}];                    % collapse lines to file
      
      if ~return_string_only
        if tt1>10, fprintf(' [write]'); end
        fwrite(out,F);
      end
      
      if tt1>10, fprintf('\n'); end
  end
  
  if ~return_string_only
    fclose(out);
    clear F
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l = slength(S)
  l=NaN;
  if isstruct(S)
    l = 0;
    if ~isempty(S) && ~isempty(fieldnames(S))
      f = fields(S);
      nf = length(f);
      len = nan(nf,1);
      for i=1:nf
        f1 = getfield(S,f{i});
        if ischar(f1), len(i) = 1;
        else len(i) = size(f1,1); end
      end
      ulen = unique(len);
      if length(ulen)==1, l = ulen;
      else
        fprintf('Warning: deprecated use of slength for structure with fields of nonuniform length\n');
        l = len(1);
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s2,ord]=sort_struct(s1,keyfield,order)
  if length(keyfield)==0, return; end
  if ~iscell(keyfield)
    keyfield = {keyfield};
  end
  if ~exist('order','var')
    order = repmat(1,length(keyfield),1);
  end
  if length(order) ~= length(keyfield)
    error('order and keyfield must have same number of elements');
  end
  
  orig_len = slength(s1);
  s2=s1;
  ord=(1:orig_len)';
  fields = fieldnames(s1);
  nf = length(fields);
  
  for k=length(keyfield):-1:1
    f = getfield(s2,keyfield{k});
    if length(f)<orig_len, error('Attempted to sort on truncated field "%s"',keyfield{k}); end
    if order(k)==1
      if isnumeric(f)
        [tmp ordi] = sortrows(f);
      else
        [tmp ordi] = sort(f);
      end
    elseif order(k)==-1
      if isnumeric(f)
        [tmp ordi] = sortrows(f,-1);
      else
        [tmp ordi] = sort(flip(f));
      end
    else
      error('Unknown order %d',order(k));
    end
    for i=1:nf
      f = getfield(s2,fields{i});
      f = f(ordi,:,:,:,:,:,:,:,:,:);
      s2 = setfield(s2,fields{i},f);
    end
    ord = ord(ordi);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = split(string, delim)
  if ischar(string), string = {string}; charflag=true; else charflag=false; end
  ns = length(string);
  R = cell(ns,1);
  for z=1:ns, if mod(z,10000)==0, fprintf('%d/%d ',z,ns); end
    
    dpos = [0 find(string{z}==delim) length(string{z})+1];
    nt = length(dpos)-1;
    R{z} = cell(nt,1);
    for t=1:nt
      R{z}{t} = string{z}(dpos(t)+1:dpos(t+1)-1);
    end
  end, if ns>=10000, fprintf('\n'); end
  if charflag, R = R{1}; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = text_to_lines(X)
  while ~isempty(X)
    if X(end)==10, X(end)=[];   % remove trailing blank lines
    else break; end
  end
  if isempty(X), L={};
  else L = split(X,char(10)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = tokens2struct(T,fieldnames)
  nt = length(T);
  nf = length(fieldnames);
  S = [];
  for f=1:nf
    F = cell(nt,1);
    for t=1:nt
      if isempty(T{t}), F{t} = '';
      else F{t} = T{t}{1}{f}; end
    end
    S = setfield(S,fieldnames{f},F);
  end
end

function a = tolines(a)
  a = split(a,char(10));
  a = a(~cellfun('isempty',a));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function verbose(str,level,varargin)

  global VERBOSE_LEVEL
  global VERBOSE_FILE

  if nargin==1
    level=1;
  end
  
  if isempty(varargin)
    
    str = char(regexprep(cellstr(str),'%','%%'));
    if ~isempty(VERBOSE_LEVEL) & (level<=VERBOSE_LEVEL)
      fprintf(1,[str repmat('\n',size(str,1),1)]');  %escape the % to prevent line from commenting
      if ~isempty(VERBOSE_FILE)
        if ~exist(VERBOSE_FILE,'file')
          fid = fopen(VERBOSE_FILE,'w');
        else
          fid = fopen(VERBOSE_FILE,'a');
        end
        
        fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
        fprintf(fid,'\n');
        fclose(fid);
      end
    end
    
  else
    
    if ~isempty(VERBOSE_LEVEL) & (level<=VERBOSE_LEVEL)
      fprintf(1,str,varargin{:});  %escape the % to prevent line from commenting
      fprintf(1,'\n')
      if ~isempty(VERBOSE_FILE)
        fid = fopen(VERBOSE_FILE,'a');
        fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
        fprintf(fid,'\n');
        fclose(fid);
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


