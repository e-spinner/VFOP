% original Value Function Optimization code altered to accept any number of
% input points and transformed into a callable function for use with
% matlab.engine by Elijah K Spinner

function result = polynomial_vfop(obj, show)

% obj consists of [firstObjectiveValue, secondObjectiveValue, rank]
% show controls if matlab plots the results or not

%Code written by Ankur Sinha
%Version: 18072011
%Note: Give non-dominated points only
%Positive epsilon value denotes that the value function correctly fitted the decision maker preferences.
       
    % IMK -- Plotting the individual points
    if show == true
        plot(obj(:,1),obj(:,2), '*r');
        hold on
    end
    
    % EKS -- Replacing square indexing horzcat with mutable size indexing
    if show == true
        for i = 1:size(obj, 1)
            label = ['P', num2str(obj(i,3))];
            text(obj(i,1), obj(i,2)+0.2, label);
        end
    end
    
    %Value function optimization. 
    %Refer: K. Deb, A. Sinha, P. Korhonen, and J.Wallenius. 
    %An interactive evolutionary multi-objective optimization method based on
    %progressively approximated value functions. IEEE Transactions on Evolutionary
    %Computation, 14(5):723739, 2010.

    % EKS -- Set optimization options
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10000);
    
    % z - covers our 
    % EKS -- passing obj to @con
    [z,fval,~,~,~] = fmincon(@fun,1.*ones(7,1),[],[],[],[],[0 0 0 0 -1000 -1000 -1000]',[1 1 1 1 1000 1000 1000]', @(x) con(x, obj), options);
    
    if fval>0
        disp('Information not correctly fitted');
    else
        disp('Information fitted correctly');
    end
    
    %    u = (.3479*x + -.1772 + .6521*y + -1.772))*(0.6521*y + -1.772 + .3479*x + -1.772);
    
    % z = [          -4.688e-04];
    
    OptimizationParameters = z; %First six are value function parameters and the last one is epsilon
    
    for i=1:size(obj,1)
        u(i)= utilfunc(obj(i,:),z);
    end
    
    for i=1:size(obj,1)
        fh = @(x,y) (z(1).*x + z(5) + z(2).*y + z(5)).*(z(3).*y + z(6) + z(4).*x + z(6)) - u(i);
    
        if show == true
            ezplot(fh, [min(obj(:,1))-1, max(obj(:,1))+1]);
        end
    end
    
    Utilities = u;
    
    if show == true
        title('');
        axis([min(obj(:,1)) max(obj(:,1)) min(obj(:,2)) max(obj(:,2))])
        xlabel('f_1')
        ylabel('f_2')
    end

    result = struct('OptimizationParameters', OptimizationParameters, 'Utilities', Utilities);
end

function f=fun(x)
    f = -x(end);
end
% Modified version 
% EKS -- Removed global obj, passing obj to constraint calculations
function [c, ceq]=con(x, obj)

    delta = 0.1;
    for i=1:size(obj,1)-1
        if obj(i,3) ~= obj(i+1,3)
            c(i) = -((utilfunc(obj(i,:),x))-(utilfunc(obj(i+1,:),x))) + x(end);
        else
            c(i) = abs((utilfunc(obj(i,:),x))-(utilfunc(obj(i+1,:),x))) - delta*x(end);
        end
    end

    %Polynomial Value Function Constraints
    % Making each term in the equation individually greater than zero %
    % IMK -- So for each P, we check each S term and make sure it's non negative
    %          We have 5 values of P (obj), and M = 2, which is each of
    %          those two lines. So really it should be a nest for-loop of P
    %          and M
    %
    % EKS -- made it iterate to double the size of obj
    for i=1:2:size(obj,1)*2
        % IMK -- This is the first term of the product in the objective
        %           function. The (i+1)/2 is just getting over the i value
        %           throughout the iteration 
        c(i+4)= -(x(1)*obj((i+1)/2,1) + x(5) + x(2)*obj((i+1)/2,2) + x(5));
        c(i+5)= -(x(3)*obj((i+1)/2,2) + x(6) + x(4)*obj((i+1)/2,1) + x(6));
    end
    % Equality constraints (Note: Minor typographic error in the paper, correct here.)

    % IMK -- one for each dimension. Makes sure terms add up to 1
    ceq(1) = -(x(1) + x(2) - 1);
    ceq(2) = -(x(3) + x(4) - 1);
end


% IMK -- Note: Utility function is the same as a value function 
% IMK -- I wrote this utility function based on eq 6 
function u=utilfunc(obj,x)
    % x(1) k11
    % x(2) k12	
    % x(3) k21	
    % x(4) k22	
    % x(5) k13	
    % x(6) k23
    u = (x(1)*obj(1) + x(5) + x(2)*obj(2) + x(5))*(x(3)*obj(2) + x(6) + x(4)*obj(1) + x(6));

end
