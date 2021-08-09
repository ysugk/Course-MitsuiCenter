function [BestX,BestLoss] = anneal_1(loss, parent, jumpsizes,options)

% ANNEAL  Minimizes a function with the method of simulated annealing
% (Kirkpatrick et al., 1983)
%
%  ANNEAL takes three input parameters, in this order:
%
%  LOSS is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single value.
%
%  PARENT is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  LT: jumpsizes(i,1) is the standard deviation of a typical jump for
%  parameter i
%
%  OPTIONS is a structure with settings for the simulated annealing. If no
%  OPTIONS structure is provided, ANNEAL uses a default structure. OPTIONS
%  can contain any or all of the following fields (missing fields are
%  filled with default values):
%
%       Verbosity: Controls output to the screen.
%                  0 suppresses all output
%                  1 gives final report only [default]
%                  2 gives temperature changes and final report 
%       Generator: Generates a new solution from an old one.
%                  Any function handle that takes a solution as input and
%                  gives a valid solution (i.e. some point in the solution
%                  space) as output.
%                  The default function generates a row vector which
%                  slightly differs from the input vector in one element:
%                  @(x) (x+(randperm(length(x))==length(x))*randn/100)
% LT: It adds a normal(0,0.01) random variable to one element of X
% LT: This one adds N(0,0.05*guess) shock to the chosen element, where guess is the initial guess value of chosen element
% LT: 		    Initialguess=parent;
% LT: 		    @(x) (x+(randperm(length(x))==length(x)).*jumpsizes*randn)
%                  Other examples of possible solution generators:
%                  @(x) (rand(3,1)): Picks a random point in the unit cube
%                  @(x) (ceil([9 5].*rand(2,1))): Picks a point in a 9-by-5
%                                                 discrete grid
%        InitTemp: The initial temperature, can be any positive number.
%                  Default is 1. 
%
%        StopTemp: Temperature at which to stop, can be any positive number
%                  smaller than InitTemp. 
%                  Default is 1e-8.
%         StopVal: Value at which to stop immediately, can be any output of
%                  LOSS that is sufficiently low for you.
%                  Default is -Inf.
%       CoolSched: Generates a new temperature from the previous one.
%                  Any function handle that takes a scalar as input and
%                  returns a smaller but positive scalar as output. 
%                  Default is @(T) (.8*T)
%      MaxConsRej: Maximum number of consecutive rejections, can be any
%                  positive number.
%                  Default is 1000. (LT: 100; after this many failed jumps, you've found solution)
%        MaxTries: Maximum number of tries within one temperature, can be
%                  any positive number. 
%                  Default is 300.
%      MaxSuccess: Maximum number of successes within one temperature, can
%                  be any positive number.
%                  Default is 20.  
%
%
%  Usage:
%     [MINIMUM,FVAL] = ANNEAL(LOSS,NEWSOL,[OPTIONS]);
%          MINIMUM is the solution which generated the smallest encountered
%          value when input into LOSS.
%          FVAL is the value of the LOSS function evaluated at MINIMUM.
%     OPTIONS = ANNEAL();
%          OPTIONS is the default options structure.
%
%
%  Example:
%     The so-called "six-hump camelback" function has several local minima
%     in the range -3<=x<=3 and -2<=y<=2. It has two global minima, namely
%     f(-0.0898,0.7126) = f(0.0898,-0.7126) = -1.0316. We can define and
%     minimise it as follows:
%          camel = @(x,y)(4-2.1*x.^2+x.^4/3).*x.^2+x.*y+4*(y.^2-1).*y.^2;
%          loss = @(p)camel(p(1),p(2));
%          [x f] = ANNEAL(loss,[0 0])
%     We get output:
%               Initial temperature:     	1
%               Final temperature:       	3.21388e-007
%               Consecutive rejections:  	1027
%               Number of function calls:	6220
%               Total final loss:        	-1.03163
%               x =
%                  -0.0899    0.7127
%               f =
%                  -1.0316
%     Which reasonably approximates the analytical global minimum (note
%     that due to randomness, your results will likely not be exactly the
%     same).

%  Reference:
%    Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
%    Simulated Annealing. _Science, 220_, 671-680.

%   joachim.vandekerckhove@psy.kuleuven.be
%   $Revision: v5 $  $Date: 2006/04/26 12:54:04 $

% Default annealing options, never used
def = struct(...
        'CoolSched',@(T) (.85*T),...
        'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
        'InitTemp',1,...
        'MaxConsRej',100,...
        'MaxSuccess',20,...
        'MaxTries',300,...
        'StopTemp',1e-8,...
        'StopVal',-Inf,...
        'Verbosity',1);

% Check input
if ~nargin %user wants default options, give it and stop
    minimum = def;
    return
elseif nargin<2, %user gave only objective function, throw error
    error('MATLAB:anneal:noParent','You need to input a first guess.');
elseif nargin<3, %user gave no options structure, use default
    options=def;
else %user gave all input, check if options structure is complete
    if ~isstruct(options)
        error('MATLAB:anneal:badOptions',...
            'Input argument ''options'' is not a structure')
    end
    fs = {'CoolSched','Generator','InitTemp','MaxConsRej',...
        'MaxSuccess','MaxTries','StopTemp','StopVal','Verbosity'};
    for nm=1:length(fs)
        if ~isfield(options,fs{nm}), options.(fs{nm}) = def.(fs{nm}); end
    end
end

% main settings
newsol = options.Generator;      % neighborhood space function
Tinit = options.InitTemp;        % initial temp
minT = options.StopTemp;         % stopping temp
cool = options.CoolSched;        % annealing schedule
minF = options.StopVal;
max_consec_rejections = options.MaxConsRej;
max_try = options.MaxTries;
max_success = options.MaxSuccess;
report = options.Verbosity;
 
k=1;

% counters etc
itry = 0;
Nbads=0; 
success = 0;
finished = 0;
consec = 0;
T = Tinit;
initenergy = loss(parent);
oldenergy = initenergy;
total = 0;
if report==2, fprintf(1,'\n  T = %7.5f, loss = %10.5f\n',T,oldenergy); end

BestX=parent; % our best X so far
BestLoss=initenergy; % best loss so far


while ~finished;
    itry = itry+1; % just an iteration counter
    current = parent; 
    
    % % Stop / decrement T criteria
    % success = number of times you've made a jump, either optimal or suboptimal
    % success >=max_success is like "too much jumping;" cool it if too much jumping
    % itry>=max_try just says "we've reached maximum trials for this temp" 
    if itry >= max_try || success >= max_success; % either finished, or cool it
        % T<minT says we've cooled enough, so give up
        % max_consec_rejections means we cannot find a better place to
        % jump, so either we've won or we've cooled enough
        if T < minT || consec >= max_consec_rejections;
            finished = 1;
            total = total + itry;
            break;
        else
        % cool when you reach maximum tries for this temp, or you're jumping too much,
        % and you haven't given up or won yet
            T = cool(T);  % decrease T according to cooling schedule                  
            if report==2, % output
                disp(['T=',num2str(T),...
                ', uprate=',num2str(Nbads/itry),...
                ', downrate=',num2str((success-Nbads)/itry),...
                ', itry=',num2str(itry),...
                ', success=',num2str(success)]);
            
                disp(['bestenergy=',num2str(BestLoss),...
                      ', Best parameters=',num2str(BestX')]);
                disp(['current energy=',num2str(oldenergy),...
                      ', current parameters=',num2str(current')]);           
                disp(' ');
            end
            total = total + itry;
            itry = 1;
            success = 1;
            Nbads=0;
        end
    end
    
    newparam = newsol(current); % new trial value
    newenergy = loss(newparam); % objective function at new trial value
    
    if (newenergy < BestLoss)
    	BestLoss=newenergy;
    	BestX=newparam;
    end
    
    if (newenergy < minF), %  We've found a winner! Stop now.
        parent = newparam; 
        oldenergy = newenergy;
        break
    end
    
    % always accept a move if it lowers the energy
    if (oldenergy-newenergy > 1e-6)  % This one has better energy, so move to it
        parent = newparam;
        oldenergy = newenergy;
        success = success+1;
        consec = 0; % we've broken our last non-better streak
    else
        if (rand < exp( (oldenergy-newenergy)/(k*T) ));  % go there anyway if we get bumped
        % more likely to bump there when newenergy is lower and when T is large
            parent = newparam;
            oldenergy = newenergy;
            success = success+1;
            Nbads=Nbads+1; % number of bad jumps
        else % stay put
            consec = consec+1; % our non-bump streak continues
        end
    end
end

%minimum = parent;
%fval = oldenergy;

if report;
    fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
    fprintf(1, '  Final temperature:       \t%g\n', T);
    fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
    fprintf(1, '  Number of function calls:\t%i\n', total);
    fprintf(1, '  Total final loss:        \t%g\n', BestLoss);
end
