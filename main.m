clear all; close all; clc;

%Setup Properties
Setup.ps = 13.3e-6;         %Pixelsize at DMD
Setup.DMDX = 256;          % number of pixels on the DMD along one axis
Setup.DMDY = 256;          % number of pixels on the DMD along the other axis
Setup.laserradius = 11e-3; %Radius of the gaussian laser spot on DMD
Setup.f = 0.1;              % focal length, f, of the telescope lens at an f-f distance between the DMD and the image plane
Setup.lambda = 635e-9;      % Wavelength of the laser light source.
N = 50;                     %Number of iterations for the optimization
P = 5;                     %number of coherent patterns being time averaged 

for z = [0.1]
    
    % generate grid
    UX = Setup.ps*(1:Setup.DMDX);
    UX = UX-mean(UX);
    UY = Setup.ps*(1:Setup.DMDY);
    UY = UY-mean(UY);
    [XX,YY] = ndgrid(UX,UY);   
    
    %laser_amplitude = ones(Setup.DMDX,Setup.DMDY); % amplitud del laser de entrada a cada espejo constante
    laser_amplitude = exp(-((XX.^2+YY.^2)/Setup.laserradius^2));
    FieldB = (laser_amplitude).*exp(1i*2*pi*rand(Setup.DMDX,Setup.DMDY));


    [NX,NY] = size(FieldB);
    psx = abs(Setup.f*Setup.lambda/(NX*Setup.ps));
    psy = abs(Setup.f*Setup.lambda/(NY*Setup.ps));

    [FieldA,~,~] = function_lens(FieldB,psx,psy,Setup.f,Setup.lambda);

    % new real grid
    %URX, URY is the real space axis limit, note the non square pixel size
    URX = psx*(1:Setup.DMDX); URX = URX-mean(URX);
    URY = psy*(1:Setup.DMDY); URY = URY-mean(URY);
    [RXX,RYY] = ndgrid(URX,URY);


    targetIntensity = mat2gray(imread('Target_Images/triangle.png'));
    target_amplitude = sqrt(targetIntensity); % saca la media cuadratica de los pixeles
    target_amplitude = target_amplitude/max(target_amplitude(:)); 
    % resize to DMD size
    target_amplitude = imresize(target_amplitude', [Setup.DMDX Setup.DMDY]);
    target_amplitudeXL = [target_amplitude rot90(target_amplitude,2)];
    % target image
    ta = double(imbinarize(imresize(target_amplitudeXL,[Setup.DMDX Setup.DMDY])));
  
    % normalize
    b = sum(ta.^2,"all");
    ta = ta/sqrt(b);
    
    % pattern initialization
    Amplitude = abs(FieldA);
    DMDPatterns = double(Amplitude>(median(Amplitude)));
    % image reconstruction with patern
    [VolumeImageAV] = function_Rendering(Setup, z,DMDPatterns);
    
    % train process
    for i = 1:N
        %Display
        f = figure(1);
        subplot(2,3,1)
        imagesc(URX,URY,ta(:,:,1)); colormap gray; axis off;
        title('Target Amplitude 1')
        subplot(2,3,2)
        imagesc(URX,URY,VolumeImageAV(:,:,1)); colormap gray; axis off;
        title('DCGH Rendering 1')


        subplot(2,3,3)
        imagesc(URX,URY,DMDPatterns); colormap gray; axis off;
        title('Pattern')

        subplot(2,3,4)
        [accuracy] = function_accuracy(VolumeImageAV,ta.^2);
        scatter(i,accuracy); hold on;
        title('Accuracy per iteration')
        ylabel('Accuracy [A.U.]')
        xlabel('Iteration [A.U.]')



        Ax_stats= subplot(2,3,5);
        set(Ax_stats,'YTickLabel',{},'XTickLabel',{},'XTick',[],'YTick',[]);
        Ax_stats.XAxis.Visible = 'off';
        Ax_stats.YAxis.Visible = 'off';
        delete(findall(gcf,'Tag','stream'));


        estimated_normalized =  rescale(VolumeImageAV);
        ta_normalized = rescale(ta);
        
        [x,d] = imhist(estimated_normalized)
        [x_,d_] = imhist(ta_normalized)

        str = {
            ['Accuracy = ', num2str(function_accuracy(VolumeImageAV,ta.^2),'%.3f')],
            ['PSNR = ', num2str(psnr(estimated_normalized(1:255,35:100),ta_normalized(1:255,35:100).^2),'%.3f')],
            ['SSIM = ', num2str(ssim(estimated_normalized(1:255,35:100),ta_normalized(1:255,35:100).^2),'%.3f')],
            ['chi-square distance =', num2str(sc_pdist2(x,x_,'cs'))]
        };

        annotation('textbox', [0.45, 0.3, 0.1, 0.1], 'String', str,'Tag','stream')
    
        %Use target amplitude to update DMD patterns
        %target amplitude = imagen objetivo
        DMDPatterns = function_GS_Z3D_inc(Setup, DMDPatterns,z,ta);
    
        %Render images from DMD patterns
        [VolumeImageAV] = function_Rendering(Setup, z,DMDPatterns);

    drawnow
    end
    
    for k = 1:P
        D = padarray(DMDPatterns,[256 384],0,'both');
        %D = DMDPatterns;
        imwrite(D,['output/' num2str(z) '_DMDPattern' '.bmp']);

    
        maxd = double(max(VolumeImageAV(:)));
        mind = double(min(VolumeImageAV(:)));
        estimate = uint8((double(VolumeImageAV) - mind)./(maxd-mind) .* 256);
        imwrite(estimate,['output/' num2str(z) '_estimate' '.png']);
    end

end
