function[frames] = animation_profile(x_array,y_array,T_struct)
figure(1)
%function that makes the animation of the plot
%first open a figure and label it 
h = mesh(x_array,y_array,T_struct(1).Temp) ;
title('Temperature Profile Against CPU Heat Generation')
xlabel('X Position (m)')
ylabel('Y Position (m)')
zlabel('T in K')
axis([0 0.01 0 0.06 293 350]) ; 

%set the number of frames equal to the amount of saved temperature profiles

Nframes = size(T_struct) ;
%allocate the struct to save the frames 
init_getframe = struct('cdata',[],'colormap',[]) ;
frames = repmat(init_getframe, Nframes(2),1) ;

%record each frame, save to a different location in the struct
for j =1:1:Nframes(2)
    set(h,'XData',x_array) ;
    set(h,'YData',y_array) ;
    set(h,'ZData',T_struct(j).Temp) ;
    drawnow
    frames(j) = getframe(gcf) ; 
    
end

%open an avi file and write to it each frame as previouslt recorded 

vidObj = VideoWriter('projectTemperatureprofile.avi');
open(vidObj);
for i=1:1:Nframes(2)
    
    writeVideo(vidObj,frames(i)) ; 
    
end
close(vidObj) ;





end