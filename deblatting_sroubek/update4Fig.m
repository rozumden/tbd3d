function update4Fig(fh,s1,s2,s3,s4)

if isempty(fh)
    return;
end
figure(fh);
if exist('s1','var') && ~isempty(s1)
    if ~isempty(s1{2})
        subplot(224);
        dispBWIm(s1{2});
    end
    if ~isempty(s1{1})
        title(s1{1});
    end
end
if exist('s2','var') && ~isempty(s2)
    if ~isempty(s2{2})
        subplot(221); dispBWIm(s2{2});
    end
    if ~isempty(s2{1})
        title(s2{1});
    end
end
if exist('s3','var') && ~isempty(s3)
    if ~isempty(s3{2})
        subplot(222); dispBWIm(s3{2});
    end
    if ~isempty(s3{1})
        title(s3{1});
    end
end 
if exist('s4','var') && ~isempty(s4)
    if ~isempty(s4{2})
        subplot(223); dispBWIm(s4{2});
    end
    if ~isempty(s4{1})
        title(s4{1});
    end
end 
drawnow;

