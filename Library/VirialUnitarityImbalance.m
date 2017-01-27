function [ muTilde1, TTilde1, impurity1, BetaMu1_vec, BetaMu2_vec, Z1_vec, Z2_vec ] = ...
    VirialUnitarityImbalance( varargin )
%VirialUnitarityImbalance This function computes the virial expansion of the EOS of 
% the imbalanced unitary fermi gas as a function of fugacity
% z=exp(\beta*\mu) of each spin component.
% (chemical pot. over Fermi energy, temperature over Fermi temperature).
% 
% Optional name value input parameters:
% "LogPoints: specifies the amount of data points on a log scale. If nothing
% is given as an input the default value is 10000.
% "Order" selects the order of the virial expansion (default value 3).
% The 'ContactOrder' can be selected between 2 and 3 (default
% value 3).

% muTilde1=mu1/E_F
% 

%addpath('/Users/Julian/Documents/MIT/MatlabPrograms/MatlabFunctionsJulian')

% Default parameters
defaultLogPoints = 5000;
defaultOrder = 3;

% create input parser
p = inputParser;

addParameter(p,'LogPoints',defaultLogPoints,@isnumeric);
addParameter(p,'Order',defaultOrder,@isnumeric)

parse(p,varargin{:});

LogPoints = p.Results.LogPoints;
order = p.Results.Order;

% Coefficient for the imbalanced virial expansion
Db2ref = 1/sqrt(2); % PHYSICAL REVIEW A 82, 023619 (2010)
Db3ref = -0.35501298; % PHYSICAL REVIEW A 82, 043626 (2010)

if order==1
    Db2=0;
    Db3=0;
elseif order==2
    Db2=Db2ref;
    Db3=0;
elseif order==3
    Db2=Db2ref;
    Db3=Db3ref;
end

%% Generate EOS data
BetaMu2_start = -20;
BetaMu2_stop = -0.01;
BetaMu1_start = -20;
BetaMu1_stop = -0.01;
BetaMu1_vec = linspace(BetaMu1_start, BetaMu1_stop, LogPoints);
Z1_vec = exp(BetaMu1_vec);
BetaMu2_vec = linspace(BetaMu2_start, BetaMu2_stop, LogPoints);
Z2_vec = exp(BetaMu2_vec);
[Z1matrix,Z2matrix]= meshgrid(Z1_vec,Z2_vec);

CommonMatrix1 = -PolyLogFrac(3/2,-Z1matrix) + 2 *Db2* Z1matrix.*Z2matrix + Db3 * (2 * Z1matrix.^2 .* Z2matrix + Z1matrix .* Z2matrix.^2);
CommonMatrix2 = -PolyLogFrac(3/2,-Z2matrix) + 2 *Db2* Z2matrix.*Z1matrix + Db3 * (2 * Z2matrix.^2 .* Z1matrix + Z2matrix .* Z1matrix.^2);


TTilde1 = 4*pi ./ (6* pi^2 * CommonMatrix1 ).^(2/3);

muTilde1 = 4*pi * log(Z1matrix) ./ (6* pi^2 * CommonMatrix1).^(2/3);

impurity1 = CommonMatrix1./CommonMatrix2;

h = surf(impurity1,TTilde1,muTilde1);
set(h,'edgecolor','none');
xlabel('n_u/n_d');
ylabel('TTilde1');
zlabel('muTilde1');
axis([0 1 0 10 -30 1])

end

