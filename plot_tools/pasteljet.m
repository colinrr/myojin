function cols = pasteljet(ncol,cpos,col0)

if nargin<1
    ncol = [];
end
if nargin<2
    cpos = [];
end
if nargin<3
    col0 = []; 
end
if isempty(ncol)
    ncol = 100;
end
if isempty(cpos)
    cpos = [0 .25 .5 .75 1];
end
if isempty(col0)
    col0 = {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f','#000000'};
end



    crgb = zeros(length(col0),3);
    for ii=1:length(col0)
        crgb(ii,:) = hex2rgb(col0{ii});
    end

    cols = zeros(ncol,3);
    ci = round(cpos.*ncol);
    for ii=1:(length(ci)-1)
        ns = diff(ci(ii:ii+1));
        cols(ci(ii)+1:ci(ii+1),1) = linspace(crgb(ii,1),crgb(ii+1,1),ns)';
        cols(ci(ii)+1:ci(ii+1),2) = linspace(crgb(ii,2),crgb(ii+1,2),ns)';
        cols(ci(ii)+1:ci(ii+1),3) = linspace(crgb(ii,3),crgb(ii+1,3),ns)';
    end
    cols = flipud(cols);
end

function rgb = hex2rgb(hexString)
	if size(hexString,2) ~= 7
		error('Not a color! %s', hexString);
	else
		r = double(hex2dec(hexString(2:3)))/255;
		g = double(hex2dec(hexString(4:5)))/255;
		b = double(hex2dec(hexString(6:7)))/255;
		rgb = [r, g, b];
	end
end
