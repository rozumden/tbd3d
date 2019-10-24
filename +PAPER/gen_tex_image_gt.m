function [srctable] = gen_tex_image_gt(gt_coeffs, img, fname, PAR, szs, ind)
scl = 1;

[~,seqname,~] = fileparts(fname);
C = strsplit(seqname,'_');
name = C{end};

scl_table = '{!}{0.2\\textwidth}';
prename = 'imgs/thumbnails/tbd/';

direc = '~/projects/vis/image';
mkdir(direc);
endtype='arrow';

src = '\\documentclass{article} \n ';

src = [src '\\usepackage[english]{babel} \n \\usepackage{graphicx} \n'];
src = [src '\\usepackage{multirow} \n'];
src = [src '\\usepackage{nopageno} \n'];
src = [src '\\usepackage{tkz-fct} \n'];
src = [src '\\usepackage{tikz,pgfplots} \n'];
src = [src '\\usetikzlibrary{arrows.meta} \n'];

src = [src '\\usepackage[margin=0.5in]{geometry} \n'];
src = [src '\n '];


src = [src '\\begin{document} \n \\noindent '];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srctable = '';

imgname = [name '.png'];
imwrite(img, fullfile(direc, imgname),'Mode','lossless');
[h,w,~] = size(img);
srctable = [srctable '\\begin{tikzpicture} \n'];
srctable = [srctable '\\begin{axis}[y dir=reverse, \n' ...
            ' xmin=1,xmax=' int2str(w) ', \n' ...
            ' ymin=1,ymax=' int2str(h) ', \n' ...
            ' xticklabels = \\empty, yticklabels = \\empty, \n' ...
            ' grid=none, axis equal image] \n'];

srctable = [srctable '\\addplot graphics[xmin=1,xmax=' int2str(w) ',ymin=1,ymax=' int2str(h) '] {' prename imgname '}; \n'];
addplot0 = '},{';
addplotE = '});  \n';

arrowstyle_default = ['>={Latex[length=1.5mm,width=0.5mm,angle''=25,open,round]},'];
lw_min_hard = 0.05;
lw_min = 0.1;
lw_max = 2;
lw_max_hard = 2.5;

parts = numel(PAR(1).R);
rr = [PAR.R]; rr = [rr(:)];
rmax = max(rr); rmin = min(rr);

dr = (rmax - rmin);
if dr / rmin < 0.1
    rmin = rmin *0.9;;
    rmax = rmax *1.1;
end

if iscell(gt_coeffs) %% Oracle
    for k = 1:numel(gt_coeffs)
        for kk = 1:numel(gt_coeffs{k})
            cf = gt_coeffs{k}{kk};
            if ~exist('szs','var') 
                rnow = PAR(k).R(kk);
            else
                xq = linspace(k, k+1, parts+1);
                rnow = interp1(ind,szs,xq(kk),'linear', 'extrap')';
            end
            alp = (rnow - rmin)./(rmax-rmin);
            lw = lw_min + alp*lw_max;
            addplotB = ['\\addplot [' arrowstyle_default ',domain=0:1,samples=2,line width=' num2str(lw) 'pt,color=yellow]({'];
            eq1 = [num2str(cf(1,1)) ' + ' num2str(cf(1,2)) '*x'];
            eq2 = [num2str(cf(2,1)) ' + ' num2str(cf(2,2)) '*x'];
            srctable = [srctable addplotB eq1 addplot0 eq2 addplotE];
        end
    end
else %% Curves
    lw_prev = lw_min;
    l_pnt = [];
    curves = gt_coeffs;
    for ki = 1:numel(curves)
        crv = curves(ki);
        % if strcmp(crv.type,'prediction')
        %     continue;
        % end
        clr = 'yellow';

        coeff = crv.coeff;
        for ci = 1:numel(coeff)
            cf = coeff{ci};
            max_power = (numel(cf) / 2) - 1;
            cf = cf*scl;
            arrow = '';
            if ci == numel(coeff) && ki == numel(curves)
                arrow = '->,';
            end
            arrowstyle = arrowstyle_default;
            if strcmp(crv.type,'bounce') || strcmp(crv.type,'connect')
                f_pnt = evaluate_coeff(cf, 0);
                if ~isempty(l_pnt) && ~all(f_pnt == l_pnt)
                    cf_ad = [l_pnt' f_pnt'];
                    cf_ad(:,2) = cf_ad(:,2) - cf_ad(:,1);
                    addplotB = ['\\addplot [' arrow arrowstyle ',domain=0:1,samples=2,line width=' num2str(lw_prev) ',color=' clr ']({'];
                    eq1 = [num2str(cf_ad(1,1)) ' + ' num2str(cf_ad(1,2)) '*x'];
                    eq2 = [num2str(cf_ad(2,1)) ' + ' num2str(cf_ad(2,2)) '*x'];
                    srctable = [srctable addplotB eq1 addplot0 eq2 addplotE];
                end

                ds = '0'; de = '1';
                ns = '2';
                lw = lw_prev;

                addplotB = ['\\addplot [' arrow arrowstyle ',domain=' ds ':' de ',samples=' ns ',line width=' num2str(lw) ',color=' clr ']({'];
                eq1 = [num2str(cf(1,1)) ' + ' num2str(cf(1,2)) '*x'];
                eq2 = [num2str(cf(2,1)) ' + ' num2str(cf(2,2)) '*x'];

                for powi = 2:max_power
                    eq1 = [eq1 ' + ' num2str(cf(1,powi+1)) '*x^' int2str(powi)];
                    eq2 = [eq2 ' + ' num2str(cf(2,powi+1)) '*x^' int2str(powi)];
                end

                srctable = [srctable addplotB eq1 addplot0 eq2 addplotE];
                
                l_pnt = evaluate_coeff(cf, 1);
            else
                f_pnt = evaluate_coeff(cf, crv.fit_iv(1));
                if ~isempty(l_pnt) && ~all(f_pnt == l_pnt)
                    cf_ad = [l_pnt' f_pnt'];
                    cf_ad(:,2) = cf_ad(:,2) - cf_ad(:,1);
                    addplotB = ['\\addplot [' arrow arrowstyle ',domain=0:1,samples=2,line width=' num2str(lw_prev) ',color=' clr ']({'];
                    eq1 = [num2str(cf_ad(1,1)) ' + ' num2str(cf_ad(1,2)) '*x'];
                    eq2 = [num2str(cf_ad(2,1)) ' + ' num2str(cf_ad(2,2)) '*x'];
                    srctable = [srctable addplotB eq1 addplot0 eq2 addplotE];
                end

                ivs = [crv.fit_iv(1):0.2:(crv.fit_iv(2)+1)];
                for kki = 2:numel(ivs)
                    kk1 = ivs(kki-1);
                    kk2 = ivs(kki);
                    if strcmp(crv.type,'joint') || strcmp(crv.type,'prediction')
                        ds = num2str(kk1);
                        de = num2str(kk2);
                        ns = int2str(ceil(10*(kk2 - kk1)));
                    end
                    kk = (kk2 + kk1)/2;
                    rnow = interp1(ind,szs,kk,'linear', 'extrap')';

                    alp = (rnow - rmin)./(rmax-rmin);
                    lw = lw_min + alp*lw_max;
                    if lw < lw_min_hard, lw = lw_min_hard; end
                    if lw > lw_max_hard, lw = lw_max_hard; end

                    lw_prev = lw;

                    addplotB = ['\\addplot [' arrow arrowstyle ',domain=' ds ':' de ',samples=' ns ',line width=' num2str(lw) ',color=' clr ']({'];
                    eq1 = [num2str(cf(1,1)) ' + ' num2str(cf(1,2)) '*x'];
                    eq2 = [num2str(cf(2,1)) ' + ' num2str(cf(2,2)) '*x'];

                    for powi = 2:max_power
                        eq1 = [eq1 ' + ' num2str(cf(1,powi+1)) '*x^' int2str(powi)];
                        eq2 = [eq2 ' + ' num2str(cf(2,powi+1)) '*x^' int2str(powi)];
                    end
                    srctable = [srctable addplotB eq1 addplot0 eq2 addplotE];
                end
                l_pnt = evaluate_coeff(cf, crv.fit_iv(2)+1);
            end
        end   
    end
end

srctable = [srctable '\\end{axis} \n'];
srctable = [srctable '\\end{tikzpicture} \n'];
srctable = [srctable '\n  \\noindent'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
src = [src srctable ' \n \\end{document} \n '];

ffile = fullfile(direc, ['image_' name '.tex']);
fid = fopen(ffile,'wt');

fprintf(fid, strrep(src,prename,''));
fclose(fid);

ffile = fullfile(direc, [name '.tex']);
fid = fopen(ffile,'wt');
srctable0 = ['\\resizebox ' scl_table ' {' srctable '}'];
fprintf(fid, srctable0);
fclose(fid);

