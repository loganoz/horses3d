function [] = printfigure(h, name, orientation,aspectratio,varargin )
if strcmp(orientation,'portrait')==1
orient portrait
elseif strcmp(orientation,'landscape')==1
    orient landscape
end
figure(h.Number);
if aspectratio == 1
    set(gcf,'PaperPositionMode', 'manual', ...
   'PaperUnits','centimeters', ...
   'Paperposition',[1 1 28 20])
    
else
    YLim = h.CurrentAxes.YLim;
    XLim = h.CurrentAxes.XLim;
    scale_factor = 1.5;
    Position = [0, 0,1280/scale_factor,int16(1280/scale_factor*aspectratio)];
    h.Position = Position;
    h.CurrentAxes.XLim = XLim;
    h.CurrentAxes.YLim = YLim;
    set(gcf,'PaperPositionMode', 'manual', ...
       'PaperUnits','centimeters', ...
       'Paperposition',Position/45)

    set(gcf,'PaperSize',[50/scale_factor,50*aspectratio/scale_factor]);
    set(gcf,'PaperPositionMode','auto');
end
if nargin>4
    handles=varargin{1};
    print(handles,'-dpdf',name);
else
    h;
    print(name,'-dpdf');
end
    
end

