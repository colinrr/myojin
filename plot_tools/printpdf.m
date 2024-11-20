function printpdf(filename,pdir,dimensions,units,dpi,renderer)
%       pdir = printpdf(filename,pdir,dimensions,units,dpi)
%   This is a quick funtion to save pdf images for figures
%   INPUT:      filename   = just that, minus extension
%               dimensions = image size, [x y]
%               path       = save directory
%               units      = units for dimensions: 'centimeters' (default)
%                                                  'inches'
%                                                  'points' (1/72")
%                                                  'normalized'
%               dpi        = dots per inch resolution [def 350]
%
% C Rowell, 2012

% To do ---
%  - varargin setup?
%  - 1 single output path is prob best? A bit funky with extension, so w/e
%  - allow [] inputs for defaulting - DONE
%

%% Parsing
narginchk(3,6)

if nargin<3 || isempty(dimensions)
    dimensions = get(gcf,'PaperPosition');
    dimensions = dimensions(3:4);
end
if nargin<4 || isempty(units)
    units = 'centimeters';
end
if nargin<5
    dpi = 350;
end
if nargin<6 || isempty(renderer)
    renderer = '-painters';
else
    assert(ismember(renderer, {'-painters','-opengl'}),'Input "renderer" must match a matlab printing renderer (see "doc print").')
end

%%

dpi = ['-r' num2str(dpi)];

set(gcf,'paperunits',units,'paperposition',[0 0 dimensions],'paperpositionmode','manual')
set(gcf,'papersize',dimensions)
fprintf('Writing figure: %s\n',fullfile(pdir, [filename '.pdf']))
print('-dpdf',dpi,fullfile(pdir, [filename '.pdf']));

