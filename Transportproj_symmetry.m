clear 
 
 
Nf = 20 ; 
%all dimensions are in meters 
FT = 0.4/1000 ; 


G = 2.2/1000 ;
L =0.008 ; %was0.051 ;


H =(48.9/1000) / 2  ;
%H =(50/1000) / 2  ;

H_cpu = (45/1000) ; 
dx = 0.0001 ;  %added /2 delete if neccessary   
dy = 0.0001 ;
%define constants and other lengths correspondinf to fin dimensions
%Initialize constants
Qmobo = 0  ;
Tmax = 337 ; %(K)
Ta = 293 ; %K
L_fin = 100/1000 ; %mm > m
W_fin = 98/1000 ; %mm > m 
k_paste = 1.7 ; %w/mK
k_cpu = 1.3 ;%W/mk
k_Cu = 385 ; %W/mK
k_Al = 205 ; %W/mK
x_cpu = 5/1000; %mm > m
x_paste = 1/1000 ; %mm > m
x_Cu = 2/1000 ; % mm > m
x_Cu_base = 1/1000 ;
h_bar= 11.38 ;%W/M^2K
h_stag = 10 ; %W/M^2K
A = W_fin*dy;
P = (2*W_fin)+(2*FT);

 
 
%begin defining bounds of mesh
M = 1 ;
x = 0 ;
ab(M) = 0 ; 

%initialize arrays for the x and y coordinates
while x< L 
    
    x = x + dx ;
    M = M + 1 ;
    ab(M) = x ; 
end

N = 1 ;
y = 0 ;
ordinate(N) = 0 ; 
while y<H
    
    y = y + dy ;
    N = N+1 ;
    ordinate(N) = y ; 
   
end

%create a mesh
[X,Y] = meshgrid(ab,ordinate) ;

%allocate an array for the temperature 
T = zeros (N,M) ;
T(:,:) = Ta ;

%initialize the average temperature to that of ambient conditions  
T_cpu_av = Ta ;
%reset variables 
x = 0 ;
y = 0 ;
m = 1 ;
n = 1 ;

%create an array that keeps track of the positions of the fins at the boundary  
fin_gap_spacing = spacer(Y,FT,G,Nf) ; 

%allocate a struct to record temperatures and against Q

Temperature_profile = struct('Temp',0,'Q',0); % note iterate this in loop to record movie 
 
p = 1 ; %index of while loop for iteration of Q
Q(p) = 1*10^7;  % Watt/s*m^3 or J

Temperature_profile(1).Temp = T ; 
Temperature_profile(1).Q = Q(p) ; 
counter = 1 ; 
i = 1 ; 

%*******
while T_cpu_av < Tmax
    %x corresponds to column and y corresponds to row in a row by column
    %mesh
    %increment Q gen here
    
    %scan all Y coordinates 
    for n = 1:1:N%y condition, scan through all y
        %scan all x coordinates  and apply appropriate boundary condition
        %based on y postion as well 
        for m=1:1:M 
            
            if X(n,m) == 0
              %left wall 
                if (Y(n,m) == 0 )
                    
                    T(n,m) = (0.5*k_cpu*T(n+1,m) + 0.5*k_cpu*T(n,m+1) + dx*h_stag*Ta + ((Q(p)*dx^2)/4) ) / (k_cpu + h_stag*dx );
               
                elseif (Y(n,m) > 0) && (Y(n,m)< H)
                    %middle of left wall
                   
                    T(n,m) =  (  (h_stag*dx*Ta + (Q(p)*dx^2)/2) + k_cpu*T(n,m+1) + 0.5*k_cpu*T(n+1,m) + 0.5* k_cpu * T(n-1,m)  ) / (2*k_cpu + h_stag*dx)  ;
                     
                 elseif (Y(n,m) < 1.0001* Y(N,m)) && (Y(n,m) > 0.9999 * Y(N,m))
                     %Top of left wall
                  
                    T(n,m) = ( h_stag*dx*0.5 + Q(p)*dx^2*0.25 + 0.5*k_cpu*T(n,m+1) + 0.5*k_cpu*T(n-1,m) ) / ( k_cpu+ 0.5*h_stag*dx) ;
                else
                    %do nothing
                end
                
            elseif (X(n,m)> 0) && (X(n,m)< x_cpu)
                
                if (Y(n,m) ==0 )% bottom of cpu surface
                    
                    T(n,m)= (0.5*k_cpu*T(n,m-1) + 0.5*k_cpu*T(n,m+1) + k_cpu*T(n+1,m) + h_stag*dx*Ta + ((Q(p)*dx^2)/2)) / (2*k_cpu + h_stag*dx) ; 
                elseif (Y(n,m) > 0) && (Y(n,m)< Y(N,m))
                  
                    T(n,m)= 0.25* (T(n,m+1) + T(n,m-1) + T(n+1,m) + T(n-1,m) + ((Q(p)*dx^2)/k_cpu) ) ; 
                elseif (Y(n,m) < 1.0001* Y(N,m)) && (Y(n,m) > 0.9999* Y(N,m)) 
                   %top of cpu surface
                  
                     T(n,m)= 0.25* (T(n,m+1) + T(n,m-1) + 2* T(n-1,m) + ((Q(p)*dx^2)/k_cpu) ) ; 
                     
                else
                    %do nothing
                end
                
            elseif (X(n,m) <1.0001* x_cpu) && (X(n,m) > 0.9999* x_cpu)
                if (Y(n,m) == 0 )
                     
                    T(n,m)=((0.5*(k_paste)*T(n,m+1))+(0.5*(k_cpu)*T(n,m-1))+ (0.5*k_cpu*T(n+1,m)) + (0.5*k_paste*T(n+1,m)) + ((h_stag)*dx)+((dx^2/4)*Q(p)))/((k_cpu)+(k_paste)+((h_stag)*dx));
                elseif (Y(n,m) > 0) && (Y(n,m)< H)
                     
                    T(n,m)= ( ((Q(p)*dx^2/2) + k_cpu*(T(n,m-1)+0.5*T(n-1,m)+0.5*T(n+1,m))) + k_paste*(T(n,m+1)+0.5*T(n-1,m)+0.5*T(n+1,m)) ) /(2*k_cpu + 2*k_paste);
                elseif (Y(n,m) < 1.0001* Y(N,m)) && (Y(n,m) > 0.9999* Y(N,m))
                    T(n,m)= ( ((Q(p)*dx^2/2) + k_cpu*(T(n,m-1)+ T(n-1,m))) + k_paste*(T(n,m+1)+ T(n-1,m)) ) /(2*k_cpu + 2*k_paste) ;
                 
                else
                    %do nothing
                end
            elseif (X(n,m)> x_cpu)&&(X(n,m)<(x_cpu + x_paste))
                if (Y(n,m) == 0 )
                    T(n,m)= ( (h_stag*dx/k_paste)*Ta + 0.5*( 2*T(n+1,m) + T(n,m+1) + T(n,m-1) )  ) / (2 + (h_stag*dx /k_paste)) ;
                elseif (Y(n,m) > 0) && (Y(n,m)< H)
                    T(n,m) = 0.25*( T(n-1,m)+T(n+1,m)+T(n,m-1)+T(n,m+1) ) ;
                elseif (Y(n,m) < 1.0001* Y(N,m)) && (Y(n,m) > 0.9999* Y(N,m))
                    T(n,m) = 0.25*( T(n-1,m) + 2*T(n,m-1)+T(n,m+1) ) ;
                    
                else
                    %do nothing
                end
            elseif (X(n,m) < 1.0001* (x_cpu + x_paste)) && (X(n,m) > 0.9999 * (x_cpu + x_paste))
                
                if (Y(n,m) == 0 )
                    T(n,m) = ((0.5*k_Cu*T(n,m+1)) + (0.5*k_paste*T(n,m-1)) + (0.5*k_paste*T(n+1,m)) + (0.5*k_Cu*T(n+1,m))  + h_stag*dx*Ta) / (k_paste + k_Cu + h_stag*dx) ;
                elseif (Y(n,m) > 0) && (Y(n,m)< H)
                    T(n,m)=(k_paste*(T(n,m-1) + 0.5*T(n+1,m) + 0.5* T(n-1,m)) + k_Cu*(T(n,m+1) + 0.5*T(n+1,m) + 0.5* T(n-1,m))  )/((2*k_paste)+(2*k_Cu));
                elseif (Y(n,m) < 1.0001* Y(N,m)) && (Y(n,m) > 0.9999* Y(N,m))
                     T(n,m)=(k_paste*(T(n,m-1) + T(n-1,m)) + k_Cu*(T(n,m+1) + T(n-1,m))  )/((2*k_paste)+(2*k_Cu));
                       
                else
                    %do nothing
                end
            elseif (X(n,m)>(x_cpu + x_paste))&&(X(n,m)< (x_cpu + x_paste + x_Cu))
                
                if (Y(n,m) == 0 )
                    T(n,m) = ( (h_stag*dx/(k_Cu)) +0.5*(2*T(n+1,m)+ T(n,m+1) + T(n,m-1))) /(2+((h_stag*dx)/k_Cu)) ;
                elseif (Y(n,m) > 0) && (Y(n,m)< H)
                    T(n,m) = 0.25*( T(n-1,m)+T(n+1,m)+T(n,m-1)+T(n,m+1) ) ;
                elseif (Y(n,m) < 1.0001* Y(N,m)) && (Y(n,m) > 0.9999* Y(N,m))
                      T(n,m) = 0.25*( 2* T(n-1,m)+T(n,m-1)+T(n,m+1) ) ;
                        
                else
                    %do nothing
                end

            elseif(X(n,m) < 1.0001*(x_cpu + x_paste + x_Cu )) && (X(n,m) > 0.9999*(x_cpu + x_paste + x_Cu ))
               
                if Y(n,m) == 0
                    
                    T(n,m) = (0.5*k_Cu*T(n,m-1) + 0.5*k_Cu*T(n+1,m) + 0.5*h_bar*dx*Ta + sqrt(h_bar*P*A*k_Al)*Ta ) / ( k_Cu + 0.5*h_bar*dx + sqrt(h_bar*P*A*k_Al) )  ; 
                elseif (Y(n,m) < 1.0001* Y(N,m)) && (Y(n,m) > 0.9999* Y(N,m))
                  
                     T(n,m) = (k_Cu*T(n,m-1) + k_Cu*T(n-1,m)  + sqrt(h_bar*P*A*k_Al)*Ta ) / (2*k_Cu + sqrt(h_bar*P*A*k_Al) )   ; 
                elseif Y(n,m)>0 && Y(n,m)<H
                    if fin_gap_spacing(n) == 1
                       % surface with fin attached
                       
                        T(n,m) = (k_Cu*T(n,m-1) + 0.5*k_Cu*T(n-1,m) + 0.5*k_Cu*T(n+1,m) + sqrt(h_bar*P*A*k_Al)*Ta ) / (2*k_Cu + sqrt(h_bar*P*A*k_Al) )   ; 
                    elseif fin_gap_spacing(n) == 0 
                       
                        %gap equations use hbar
                      T(n,m) = (k_Cu* T(n,m-1) + 0.5*k_Cu* T(n+1,m) + 0.5*k_Cu*T(n-1,m) + h_bar*dx*Ta  ) / (2*k_Cu + h_bar*dx )  ; 
                        
                    else
                        
                    end
                    
                else
                   %nothing  
                end
            end
        end
    end
% increment loop counter 
p = p + 1 ; 
%increment generation of heat 
Q(p) = Q(p-1) +1000  ; 
%set cpu bound 
cpu_bound =  find(X(1,:)== x_cpu) ; 
%condense to an array that records only the CPU temp 
Tcpu = T(:,1:cpu_bound) ;

%account for the symmetry adiabat 
Thalf = flipud(T);
Tfull = [T ; Thalf] ; 

%********purposes of creating animation
if p == 1 
    %record the temperature profiles at different Q values and the if
    %statements account for indexin of an odd number . Must increment in
    %large steps to save memory and avoid recording too many frames 
Temperature_profile(i).Temp = Tfull ; 
Temperature_profile(i).Q = Q(p) ;
i = i + 1 ;

elseif (p > 1) && (p == 10*counter)
Temperature_profile(i).Temp = Tfull ; 
Temperature_profile(i).Q = Q(p) ;
counter = counter + 1 ; 
i = i + 1 ;
elseif rem(p,10) > 0
Temperature_profile(i).Temp = Tfull ; 
Temperature_profile(i).Q = Q(p) ;
end 
%Evaluate the average temperature as the loop condition 
T_cpu_av = mean(mean(Tcpu)) ;


%*********
end

%create a new mesh for the full profile reflected over the line of symmetry
%in the middle 
x = 0 ; 
ab(1) = 0 ; 

for i=2:1:M
  x = x + dx ; 
  ab(i) = x ; 
     
end


y = 0 ;
ordinate(1) = 0 ; 

for i=2:1:2*size(Thalf(:,1))
 
    y = y + dy ;
    ordinate(i) = y ; 
  
     
end

[X_plot,Y_plot] = meshgrid(ab,ordinate) ;
%call to the animation function I created
animation_profile(X_plot,Y_plot,Temperature_profile) 
%calculate the wattage of heat generated in the chip by multiplying by its
%volume in meters 
Qmax_rate = Q(p) *0.0489*0.0489*0.005 ; 
%print the max heat value with appropriate units 
fprintf('\nThe Maximum CPU heat generation rate is %5.2f Watts',Qmax_rate)
