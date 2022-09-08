%
%////////////////////////////////////////////////////////////////
%
%	This code evaluates the Kinetic Energy spectra
%	of a given 3D velocity field.
%
%/////////////////////////////////////////////////////////////////
%
%clc; clear all; close all;
%
%   Program options
%   ---------------
    cubeFile = '/media/magerit/horses_cases/TaylorGreen/PaperReviewComputations/RESULTS/TGV8P2Psvv01/TGV8P2Psvv01_0000001602.cube.tec';
    plotSpectra = true;
    writeSpectra = false;
    resultsName = 'spectra';
%
%   Read raw data from file
%   -----------------------
    rawData = importdata(cubeFile);
%
%   Number of points in each direction
%   ----------------------------------
    Np=round(size(rawData.data,1)^(1/3) -1);
%
%   Cube length
%   -----------
    L=max(rawData.data(:,1)) - min(rawData.data(:,1));
%
%   Minimal and Nyquist wave number
%   -------------------------------
    k0   = 2.0*pi/L;
    kNYQ = sqrt(3)*0.5*Np*k0;
%
%   Wave numbers
%   ------------
    kx=k0*(-Np/2:1:Np/2-1);
    ky=k0*(-Np/2:1:Np/2-1);
    kz=k0*(-Np/2:1:Np/2-1);
%
%   Reshape the vectors to gather cartesian velocity components
%   -----------------------------------------------------------
    ux=reshape(rawData.data(:,4),Np+1,Np+1,Np+1);
    uy=reshape(rawData.data(:,5),Np+1,Np+1,Np+1);
    uz=reshape(rawData.data(:,6),Np+1,Np+1,Np+1);
    
    ux = ux(1:Np,1:Np,1:Np);
    uy = uy(1:Np,1:Np,1:Np);
    uz = uz(1:Np,1:Np,1:Np);
%
%   Get the Kinetic Energy
%   ----------------------
    KinEn = 0.5 * sum(sum(sum(((ux.*ux + uy.*uy + uz.*uz)))))/(L^3);
%
%   Perform 3D Fourier transforms (fft)
%   -----------------------------------
    Ux=fftn(ux)/(Np^3);
    Uy=fftn(uy)/(Np^3);
    Uz=fftn(uz)/(Np^3);
%
%   Shift in wavenumber increasing order
%   ------------------------------------
    Ux=fftshift(Ux);
    Uy=fftshift(Uy);    
    Uz=fftshift(Uz);    
%
%   ******************************************************************
%   Compute the 3D Spectra: energy is computed in spherical crowns of
%   radius          (n-0.5)*k0 < k3d < (n+0.5)*k0
%   and stored in the E(n) arrays.
%   ******************************************************************
%
%   Get the 1D wavenumbers: 3D wavenumbers range up to sqrt(3)k0(Np/2)
%   Hence, we specify k0*Np to make sure we encompass all 3D wavenumbers
%   --------------------------------------------------------------------
    k_1d=k0.*(0:Np);
%
%   Initialization  
%   --------------
    E=zeros(Np+1,1);
    
    for i=1:Np ; for j=1:Np ; for m=1:Np 
%
%      Get the 3D wavenumber
%      ---------------------
       k_abs =sqrt(kx(i)^2+ky(j)^2+kz(m)^2 );
%
%      Get the crown and append the energy
%      -----------------------------------
       n = floor(k_abs/k0 + 0.5);
       E(n+1)=E(n+1)+0.5*(abs(Ux(i,j,m))^2+abs(Uy(i,j,m))^2+abs(Uz(i,j,m))^2);
           
    end ;       end ;         end 
%
%   Plot results
%   ------------
    if ( plotSpectra ) 
        h = figure(1);
        hold on
        plot(k_1d,E,'-k');
        %plot(k_1d,E/KinEn,'ok','MarkerSize',8,'MarkerFaceColor',[1,1,1]);
        e_kolmogorov=E(3)*(k_1d/k_1d(3)).^(-5/3);
        plot(k_1d,e_kolmogorov,'--m','LineWidth',1);
        plot([kNYQ,kNYQ],[1e-20,1e0],'--r');
        axis([k0,2*kNYQ,1e-12,e_kolmogorov(2)]);
%
%       Figure options
%       --------------
        h.CurrentAxes.FontSize = 20;
        h.CurrentAxes.YScale = 'log';
        h.CurrentAxes.XScale = 'log';
        h.CurrentAxes.YMinorTick = 'on';
        h.CurrentAxes.XMinorTick = 'on';
        h.CurrentAxes.LineWidth = 2;
        h.CurrentAxes.YTick = 10.^(-16:4:3);
        %h.CurrentAxes.YTickLabel = 
        %set(haxes(1),'YTick',0:0.001:1)
        xlabel('$k$');
        ylabel('$E(k)$');
        printfigure(h,resultsName,'landscape',0.3);
    end
%
%   Write results to file
%   ---------------------
    if ( writeSpectra )     
        spectra=[transpose(k_1d), E, transpose(e_kolmogorov)];
        save([resultsName,'.dat'],'spectra','-ASCII');
    end

    



