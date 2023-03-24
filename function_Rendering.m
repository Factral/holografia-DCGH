function [VolumeImage] = function_Rendering(Setup, z,DMDPatterns)
%This function renders images at distance z from DMD patterns
%Inputs:
%z = an array of size LZ corresponding to the propagation distance of the image depth planes
%DMDPatterns = a set of p DMD pattterns of size LX, by LY, by p
%Outputs:
%VolumeImage = A stack of intensity data corresponding to individually propagated oherent waves, of size 
%LX by LY  by LZ by p

    UX = Setup.ps*(1:Setup.DMDX); UX = UX-mean(UX);
    UY = Setup.ps*(1:Setup.DMDY); UY = UY-mean(UY);
    [XX,YY] = ndgrid(UX,UY);
    %
    
    % amplitud =1
    %laser_amplitude = ones(Setup.DMDX,Setup.DMDY);
    laser_amplitude = exp(-((XX.^2+YY.^2)/Setup.laserradius^2));
    %
    
    %
    mask = double(sqrt(XX.^2+YY.^2)<0.0001);  % mascara del centro de la grilla
    %
    
    VolumeImage = zeros(Setup.DMDX,Setup.DMDY,1,1); % inicializa matriz
    
    % como laseramplitude = 1
    % es innecesario este paso
    FieldA = laser_amplitude.*DMDPatterns;
    
    %go to image plane
    [FieldB,psx,psy] = function_lens(FieldA,Setup.ps,Setup.ps,-Setup.f,Setup.lambda);
    
    %normaliza
    FieldB = FieldB/sqrt((sum(abs(FieldB(:).^2))));
    
    fieldZ = function_propagate((1-mask).*FieldB,Setup.lambda,z,psy, psx);
    
    fieldZ = fieldZ/sqrt((sum(abs(fieldZ(:).^2))));
    VolumeImage(:,:,1,1) = abs(fieldZ.^2);

end