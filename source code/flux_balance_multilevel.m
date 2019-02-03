

% FLUX_BALANCE Multi-level flux-balance analysis of FBA model for multiple biomass objectives in a gut community
%    [V, FMAX, FMIN] = FLUX_BALANCE(MODELJOINT) performs a basic flux-balance
%    analysis of the FBA model MODELJOINT and returns a biomass-maximizing
%    flux distribution in the vector V.  The maximum and minimum 
%    synthetic objective possible in a biomass-maximizing flux distribution
%    is given in FMAX and FMIN, respectively.
%
%    [V, FMAX, FMIN] = FLUX_BALANCE(MODELJOINT, QUIET) performs the analysis
%    and suppresses screen output if QUIET is set to true.

function [v10max_max, f10max_max] = flux_balance_multilevel(modelJoint)

% if nargin < 2
%     quiet = false;
% end

param.tmlim  = -1;
param.msglev = 1;
param.save   = 0;

nrxn   = numel(modelJoint.rxns);
nmetab = numel(modelJoint.mets);

yt = ones(nrxn,1); %the old yt = modelJoint.present, meaning that all the reactions must be considered as active in the model

A = [ modelJoint.S; 
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1) ];
ctype = char('S' * ones(1, nmetab + nnz(~yt)));
vartype = char('C' * ones(1, nrxn));
[v, vbiomass] = glpk(modelJoint.f, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);


%min and max 2nd objective modelJoint.g
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f' ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass ];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 1));     %An array of characters containing the sense of each constraint in the constraint matrix.  Each element of the array may be one of the following values %           'F' Free (unbounded) variable (the constraint is ignored).  %           'U' Variable with upper bound ( A(i,:)*x <= b(i)).  %           'S' Fixed Variable (A(i,:)*x = b(i)).   %           'L' Variable with lower bound (A(i,:)*x >= b(i)).   %           'D' Double-bounded variable (A(i,:)*x >= -b(i) and A(i,:)*x <= b(i)).
%[v1min, fmin] = glpk(modelJoint.g, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v1max, fmax] = glpk(modelJoint.g, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

%min and max 3rd objective modelJoint.h with minimum 2nd objective modelJoint.g
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 2));     %one 'S' more than the case before because there is another constraint
% [v1min_min, fmin_min] = glpk(modelJoint.h, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v1min_max, fmin_max] = glpk(modelJoint.h, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

%min and max 3rd objective modelJoint.h with maximum 2nd objective modelJoint.g
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 2)); 
%[v1max_min, fmax_min] = glpk(modelJoint.h, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v1max_max, fmax_max] = glpk(modelJoint.h, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 4th objective modelJoint.i with minimum 3rd objective modelJoint.h
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
      % fmin_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 3));     %one 'S' more than the case before because there is another constraint
% [v2min_min, f2min_min] = glpk(modelJoint.i, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v2min_max, f2min_max] = glpk(modelJoint.i, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 4th objective modelJoint.i with maximum 3rd objective modelJoint.h
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
      fmax_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 3));     %one 'S' more than the case before because there is another constraint
%[v2max_min, f2max_min] = glpk(modelJoint.i, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v2max_max, f2max_max] = glpk(modelJoint.i, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

%min and max 5th objective modelJoint.j with minimum 4th objective modelJoint.i
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 4));     %one 'S' more than the case before because there is another constraint
% [v3min_min, f3min_min] = glpk(modelJoint.j, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v3min_max, f3min_max] = glpk(modelJoint.j, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 5th objective modelJoint.j with maximum 4th objective modelJoint.i
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
	  fmax_max;
	  f2max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 4));     %one 'S' more than the case before because there is another constraint
%[v3max_min, f3max_min] = glpk(modelJoint.j, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v3max_max, f3max_max] = glpk(modelJoint.j, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 6th objective modelJoint.k with minimum 5th objective modelJoint.j
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i';
	  % modelJoint.j'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min;
	  % f3min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 5));     %one 'S' more than the case before because there is another constraint
% [v4min_min, f4min_min] = glpk(modelJoint.k, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v4min_max, f4min_max] = glpk(modelJoint.k, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 6th objective modelJoint.k with maximum 5th objective modelJoint.j
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i';
	  modelJoint.j'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
	  fmax_max;
	  f2max_max;
	  f3max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 5));     %one 'S' more than the case before because there is another constraint
%[v4max_min, f4max_min] = glpk(modelJoint.k, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v4max_max, f4max_max] = glpk(modelJoint.k, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 7th objective modelJoint.l with minimum 6th objective modelJoint.k
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i';
	  % modelJoint.j';
	  % modelJoint.k'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min;
	  % f3min_min;
	  % f4min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 6));     %one 'S' more than the case before because there is another constraint
% [v5min_min, f5min_min] = glpk(modelJoint.l, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v5min_max, f5min_max] = glpk(modelJoint.l, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 7th objective modelJoint.l with maximum 6th objective modelJoint.k
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i';
	  modelJoint.j';
	  modelJoint.k'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
	  fmax_max;
	  f2max_max;
	  f3max_max;
	  f4max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 6));     %one 'S' more than the case before because there is another constraint
%[v5max_min, f5max_min] = glpk(modelJoint.l, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v5max_max, f5max_max] = glpk(modelJoint.l, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 8th objective modelJoint.m with minimum 7th objective modelJoint.l
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i';
	  % modelJoint.j';
	  % modelJoint.k';
	  % modelJoint.l'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min;
	  % f3min_min;
	  % f4min_min;
	  % f5min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 7));     %one 'S' more than the case before because there is another constraint
% [v6min_min, f6min_min] = glpk(modelJoint.m, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v6min_max, f6min_max] = glpk(modelJoint.m, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 8th objective modelJoint.m with maximum 7th objective modelJoint.l
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i';
	  modelJoint.j';
	  modelJoint.k';
	  modelJoint.l'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
	  fmax_max;
	  f2max_max;
	  f3max_max;
	  f4max_max;
	  f5max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 7));     %one 'S' more than the case before because there is another constraint
%[v6max_min, f6max_min] = glpk(modelJoint.m, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v6max_max, f6max_max] = glpk(modelJoint.m, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

%min and max 9th objective modelJoint.n with minimum 8th objective modelJoint.m
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i';
	  % modelJoint.j';
	  % modelJoint.k';
	  % modelJoint.l';
	  % modelJoint.m'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min;
	  % f3min_min;
	  % f4min_min;
	  % f5min_min;
	  % f6min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 8));     %one 'S' more than the case before because there is another constraint
% [v7min_min, f7min_min] = glpk(modelJoint.n, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v7min_max, f7min_max] = glpk(modelJoint.n, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 9th objective modelJoint.n with maximum 8th objective modelJoint.m
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i';
	  modelJoint.j';
	  modelJoint.k';
	  modelJoint.l';
	  modelJoint.m'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
	  fmax;
	  fmax_max;
	  f2max_max;
	  f3max_max;
	  f4max_max;
	  f5max_max;
	  f6max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 8));     %one 'S' more than the case before because there is another constraint
%[v7max_min, f7max_min] = glpk(modelJoint.n, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v7max_max, f7max_max] = glpk(modelJoint.n, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

%min and max 10th objective modelJoint.o with minimum 9th objective modelJoint.n
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i';
	  % modelJoint.j';
	  % modelJoint.k';
	  % modelJoint.l';
	  % modelJoint.m';
	  % modelJoint.n'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min;
	  % f3min_min;
	  % f4min_min;
	  % f5min_min;
	  % f6min_min;
	  % f7min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 9));     %one 'S' more than the case before because there is another constraint
% [v8min_min, f8min_min] = glpk(modelJoint.o, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v8min_max, f8min_max] = glpk(modelJoint.o, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 10th objective modelJoint.o with maximum 9th objective modelJoint.n
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i';
	  modelJoint.j';
	  modelJoint.k';
	  modelJoint.l';
	  modelJoint.m';
	  modelJoint.n'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
	  fmax_max;
	  f2max_max;
	  f3max_max;
	  f4max_max;
	  f5max_max;
	  f6max_max;
	  f7max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 9));     %one 'S' more than the case before because there is another constraint
%[v8max_min, f8max_min] = glpk(modelJoint.o, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v8max_max, f8max_max] = glpk(modelJoint.o, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 11th objective modelJoint.p with minimum 10th objective modelJoint.o
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i';
	  % modelJoint.j';
	  % modelJoint.k';
	  % modelJoint.l';
	  % modelJoint.m';
	  % modelJoint.n';
	  % modelJoint.o'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min;
	  % f3min_min;
	  % f4min_min;
	  % f5min_min;
	  % f6min_min;
	  % f7min_min;
	  % f8min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 10));     %one 'S' more than the case before because there is another constraint
% [v9min_min, f9min_min] = glpk(modelJoint.p, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v9min_max, f9min_max] = glpk(modelJoint.p, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 11th objective modelJoint.p with maximum 10th objective modelJoint.o
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i';
	  modelJoint.j';
	  modelJoint.k';
	  modelJoint.l';
	  modelJoint.m';
	  modelJoint.n';
	  modelJoint.o'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
	  fmax_max;
	  f2max_max;
	  f3max_max;
	  f4max_max;
	  f5max_max;
	  f6max_max;
	  f7max_max;
	  f8max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 10));     %one 'S' more than the case before because there is another constraint
%[v9max_min, f9max_min] = glpk(modelJoint.p, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v9max_max, f9max_max] = glpk(modelJoint.p, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 12th objective modelJoint.q with minimum 11th objective modelJoint.p
% A = [ modelJoint.S;
      % sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      % modelJoint.f';
      % modelJoint.g';
	  % modelJoint.h';
	  % modelJoint.i';
	  % modelJoint.j';
	  % modelJoint.k';
	  % modelJoint.l';
	  % modelJoint.m';
	  % modelJoint.n';
	  % modelJoint.o';
	  % modelJoint.p'];
% b = [ zeros(nmetab, 1);
      % zeros(nnz(~yt), 1);
      % vbiomass;
      % fmin;
	  % fmin_min;
	  % f2min_min;
	  % f3min_min;
	  % f4min_min;
	  % f5min_min;
	  % f6min_min;
	  % f7min_min;
	  % f8min_min;
	  % f9min_min];
% ctype = char('S' * ones(1, nmetab + nnz(~yt) + 11));     %one 'S' more than the case before because there is another constraint
% [v10min_min, f10min_min] = glpk(modelJoint.q, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
% [v10min_max, f10min_max] = glpk(modelJoint.q, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);

% min and max 12th objective modelJoint.q with maximum 11th objective modelJoint.p
A = [ modelJoint.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      modelJoint.f';
      modelJoint.g';
	  modelJoint.h';
	  modelJoint.i';
	  modelJoint.j';
	  modelJoint.k';
	  modelJoint.l';
	  modelJoint.m';
	  modelJoint.n';
	  modelJoint.o';
	  modelJoint.p'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax;
	  fmax_max;
	  f2max_max;
	  f3max_max;
	  f4max_max;
	  f5max_max;
	  f6max_max;
	  f7max_max;
	  f8max_max;
	  f9max_max];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 11));     %one 'S' more than the case before because there is another constraint
%[v10max_min, f10max_min] = glpk(modelJoint.q, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, 1);
[v10max_max, f10max_max] = glpk(modelJoint.q, A, b, modelJoint.lb, modelJoint.ub, ctype, vartype, -1);


% if ~quiet
%     fprintf('Biomass flux          %s:    %f\n',modelJoint.rxns{find(modelJoint.f==1)}, modelJoint.f' * v)
%     fprintf('1st Synthetic flux    %s:  [fmax = %.15f]\n', modelJoint.rxns{find(modelJoint.g==1)}, fmax)
%     fprintf('2nd Synthetic flux    %s:  [fmax_max = %.15f]\n', modelJoint.rxns{find(modelJoint.h==1)}, fmax_max)
% 	  fprintf('3rd Synthetic flux    %s:  [f2max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.i==1)}, f2max_max)
% 	  fprintf('4th Synthetic flux    %s:  [f3max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.j==1)}, f3max_max)
%     fprintf('5th Synthetic flux    %s:  [f4max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.k==1)}, f4max_max)
%     fprintf('6th Synthetic flux    %s:  [f5max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.l==1)}, f5max_max)
%     fprintf('7th Synthetic flux    %s:  [f6max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.m==1)}, f6max_max)
%     fprintf('8th Synthetic flux    %s:  [f7max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.n==1)}, f7max_max)
%     fprintf('9th Synthetic flux    %s:  [f8max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.o==1)}, f8max_max)
%     fprintf('10th Synthetic flux   %s:  [f9max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.p==1)}, f9max_max)
%     fprintf('11th Synthetic flux   %s:  [f10max_max = %.15f]\n', modelJoint.rxns{find(modelJoint.q==1)}, f10max_max)
% 	end