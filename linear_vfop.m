%Code written by Ankur Sinha
%Version: 18072011
%Note: Give non-dominated points only
%Positive epsilon value denotes that the value function correctly fitted the decision maker preferences.

function result = linear_vfop( obj, show )

    %Example Partial ordering
    %obj = [3.6 3.9 1;     2.5 4.1 2;   5.5 2.5 3;     0.5 5.2 3;    6.9 1.8 4];
    
    if show == true
        plot(obj(:,1),obj(:,2), '*r');
        hold on
    
        arrow = [];
        for i=1:size(obj,1)
            %arrow(i,:) = '\leftarrow';
            arrow(i,:) = 'P';
        end
        
        text(obj(:,1),obj(:,2)+0.2,[arrow(:,:) num2str(obj(:,3))]);
    end
    
    %Value function optimization. 
    %Refer: K. Deb, A. Sinha, P. Korhonen, and J.Wallenius. 
    %An interactive evolutionary multi-objective optimization method based on
    %progressively approximated value functions. IEEE Transactions on Evolutionary
    %Computation, 14(5):723–739, 2010.
    
    % EKS -- Set optimization options
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10000);
    
    % EKS -- passing obj to @con
    [z,fval,~,~,~] = fmincon(@fun,0.5.*ones(3,1),[],[],[],[],[0 0 -1000]',[1 1 1000]', @(x) con(x, obj), options);
    
    if fval>0 
        disp('Information not correctly fitted');
    else
        disp('Information fitted correctly');
    end
    
    OptimizationParameters = z; %First six are value function parameters and the last one is epsilon
    
    for i=1:size(obj,1)
        u(i)= utilfunc(obj(i,:),z);
    end
    
    if show == true
        
        for i=1:size(obj,1)
            fh = @(x,y) (x*z(1)+y*z(2)) - u(i);
            ezplot(fh, [min(obj(:,1))-1, max(obj(:,1))+1]);
        end
    
    
        title('');
        axis([min(obj(:,1)) max(obj(:,1)) min(obj(:,2)) max(obj(:,2))])
        xlabel('f_1')
        ylabel('f_2')
    end
    
    Utilities = u;
    
    result = struct('OptimizationParameters', OptimizationParameters, 'Utilities', Utilities);

end

%% Functions
% start -- Objective function that simply inverts the ranking value in the input 
function f=fun(x)

  f = -x(end);

end
% end   -- Objective function 

% start -- Constraint function
function [c, ceq]=con(x, obj)

  delta = 0.1;
  for i=1:size(obj,1)-1
      if obj(i,3) ~= obj(i+1,3)
          c(i) = -((utilfunc(obj(i,:),x))-(utilfunc(obj(i+1,:),x))) + x(end); % IMK - x(end) is epsilon 

          % Minimization version
          %c(i) = -((utilfunc(obj(i+1,:),x))-(utilfunc(obj(i,:),x))) + x(end); % IMK - x(end) is epsilon 


      else
          c(i) = abs((utilfunc(obj(i,:),x))-(utilfunc(obj(i+1,:),x))) - delta*x(end); % IMK - x(end) is epsilon 
      end
  end
  ceq(1)=x(1) + x(2) - 1;

end
% end   -- Constraint function

% start -- Utility function
function u=utilfunc(obj,x)
    u = obj(1)*x(1)+obj(2)*x(2);
end
% end   -- Utility function


