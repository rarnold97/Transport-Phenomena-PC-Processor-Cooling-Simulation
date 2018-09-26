function[Fins] = spacer(Y_values,Thickness,Gap_thickness, N_fin) %pass in the Ycorrdinates of the edge, the gapthickness, fin thickness, and number of fins 
i = 1 ; 
% function that determines where fins are 
while Y_values(i) < (Thickness + Gap_thickness)
    
    if Y_values(i) < Thickness
        %if withing one gap thickness, record a 1
        array(i) = 1 ;
        
    else
        %if in a gap record a zero
        array(i) = 0 ; 
        
    end
%increment loop counter 
i = i + 1 ; 



end

%transpose in to a row vecotr and copy array Nfin + 1 times to account for the fact that there are more fins that gaps if the top and bottom of the heat sink have a fin border 
Fins = transpose(array) ; 
Fins = repmat(array,1,(N_fin+1));

end