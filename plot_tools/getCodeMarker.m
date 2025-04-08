function [codeMarker, codeDesc, codeAlpha] = getCodeMarker(codeSub)
% codeMarkers = {'o','s','^'};
% sCodeDescriptions = ["Intrusive","Effusive","Fragmenting"];


    if codeSub==3 || codeSub==2 % P-bal and explosive frag
        codeDesc = "Fragmenting";
        codeMarker = '^';
        codeAlpha = 1.0;
    elseif codeSub==1 || codeSub==-1 % Taking all effusive as valid (see Obsidian notes)
        codeDesc = "Effusive";
        codeMarker = 's';
        codeAlpha = 0.8;
    elseif codeSub==0
        codeDesc = "Intrusive";
        codeMarker = 'o';
        codeAlpha = 0.6;
    else
        codeDesc = "Invalid";
        codeMarker = '.';
        codeAlpha = 0.25;
    end
    
%     alpha_max = 1;
%     alpha_min = 0.25;
    
end