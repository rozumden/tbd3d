function [curves, frms] = sequence_fit(frms, im, expos, power_factor, remove_single_fmod)
if nargin < 5
	remove_single_fmod = true;
end

if nargin < 4
	power_factor = 3;
end

if nargin < 3
	expos = 1;
end

max_outliers = 1;

for id = 1:numel(frms)
	if ~isempty(frms{id})
		frms{id}.instance = id;
	end
end

st_id = 1;
en_id = 0;
intervals = {};
was_tbd = true; last_fmod_id = 0;
non_pred_id = ~cellfun(@(x) ~isempty(x) && strcmp(x.caseused, 'Pred'), frms);
for id = 1:numel(frms)
	frm = frms{id};
	if isempty(frm) || strcmp(frm.caseused, 'Pred')
		if sum(non_pred_id(id:end)) > 0
			frms{id} = [];
		end
		continue;
	end

	if ~isempty(frm.Direction)
		if abs(frm.Direction) == 1
			[st1, en1] = FIT.get_sten(frm);
			frm.End = en1';
			frm.Start = st1';
		end
	else
		frm.End = frm.Centroid;
	end

	if strcmp(frm.caseused, 'TbD') || ~remove_single_fmod
		was_tbd = true;
		last_fmod_id = 0;
	elseif strcmp(frm.caseused, 'FMOd')
		if ~was_tbd
			frms{last_fmod_id} = [];
		end
		was_tbd = false;
		last_fmod_id = id;
	end

	frms0 = [frms{st_id:id-1}];
	if isempty(frms0), continue; end
	d_to_st = abs(frms0(1).End - frm.End);
	[~,dind] = max(d_to_st);
	d = arrayfun(@(x) frm.End(dind) - x.End(dind), frms0,'UniformOutput',false);
	d = reshape([d{:}], 1, numel(frms0));
	ds = sign(d); %% signs
	dm = mode(ds')'; %% most common sign
	dout = sum([ds ~= dm]')'; %% number of outliers
	%% keep last known perfect interval
	if min(dout) == 0 && (~strcmp(frm.caseused, 'FMOd')  || ~remove_single_fmod)
		en_id = id; 
	end

	if min(dout) > max_outliers || id == numel(frms)
		if en_id == 0, error('No perfect intervals!'); end
		intervals = [intervals {[st_id en_id]}];
		st_id = en_id+1;
		en_id = 0;
	end
end

if en_id ~= 0
	intervals = [intervals {[st_id numel(frms)]}];
end
for kkk = 1:numel(frms)
	if ~isempty(frms{kkk}) && isempty(frms{kkk}.Direction)
		frms{kkk} = [];
	end
end

%% Split intervals
min_curve_length = 10;
wsize = 13;
curves = [];
for ivid = 1:numel(intervals)
	iv = intervals{ivid};
	fs = frms(iv(1):iv(2));
	[curve] = FIT.sequence_fit_part(fs, im);
	dc = [[0; 0] diff(curve')']; 
	[~,dcind] = min( sum(dc' == mode(dc')) );
	c1 = curve(dcind,:);
	ddcind = zeros(size(c1));
	for kki = wsize+1:numel(c1)-wsize-1
		nbr = c1(kki-wsize:kki+wsize);
		el = nbr(wsize+1);
		dfv = nbr - el;
		dfall = [0 diff(nbr)];
		if numel(unique(sign(dfv))) == 3, continue; end %% not extreme value
		if max(abs(dfv)) <= 3, continue; end %% flat region
		%% perfectly flat on one side at least
		if max(abs(dfv(1:wsize))) < 1 || max(abs(dfv(wsize+2:end))) < 1, continue; end 
		s1 = unique(sign(dfall(1:wsize))); s1 = s1(s1 ~= 0);
		s2 = unique(sign(dfall(wsize+2:end))); s2 = s2(s2 ~= 0);
		if numel(s1) > 1 || numel(s2) > 1, continue; end %% left or right is not monotonous
		ddcind(kki) = 1;
	end
	%% remove small regions
	for kki = find(ddcind)
		[~,last_ind] = find(ddcind(1:kki-1)); 
		if isempty(last_ind), continue; end
		last_ind = last_ind(end);
		if kki - last_ind < min_curve_length
			ddcind(last_ind:kki) = 1;
		end
	end

	%% remove clusters
	dt = diff(bwdist(~ddcind));
	temp = [0 dt; dt 0];
	break_points = ismember(temp',[1 -1], 'rows')' | ismember(temp',[1 0], 'rows')';

	break_points(1) = 1;
	break_points(end) = 1;
	break_points = find(break_points);
	ix = FIT.frms2curve(fs, curve);
	for kki = 2:numel(break_points)
		crv = [];
		crv.iv = iv;
		crv.ivid = ivid;
		crv.curve = curve(:,break_points(kki-1):break_points(kki));
		ixinl = ix >= break_points(kki-1) & ix <= break_points(kki);
		full_in_curve = find(sum(ixinl) == 2); 
		if isempty(full_in_curve)
			continue;
		else
			st_fit = full_in_curve(1);
			en_fit = full_in_curve(end);
			max_power = max(2, min(6, round((en_fit-st_fit)/power_factor) ));
			if en_fit-st_fit == 0, max_power = 1; end
			coeff = FIT.lsq_fit(frms, st_fit, en_fit, max_power, expos);
		end
		[crv.Fit,~,crv.len] = myTrajRender(size2(im), coeff, [st_fit en_fit+1]);
		crv.coeff = coeff;
		crv.fit_iv = [st_fit en_fit];
		crv.type='joint';
		curves = [curves crv];

		%% predict to past
		if kki == 2 && ivid == 1
			if curves(1).fit_iv(1) > 1 
				crv = [];
				crv.iv = curves(1).iv;
				crv.ivid = curves(1).ivid;
				crv.coeff = curves(1).coeff;
				crv.fit_iv = [1 curves(1).fit_iv(1)-1];
				crv.curve = curves(1).curve;
				[crv.Fit,~,crv.len] = myTrajRender(size2(im), crv.coeff, [1 crv.fit_iv(2)+1]);
				crv.type = 'prediction';
				curves = [crv curves];
			end
		end
		%% start reconstructing holes
		if kki == 2 && ivid > 1 && numel(curves) > 1 %% connect previous curves 
			crv0 = curves(end-1);
			in0 = crv0.fit_iv(2);
			in1 = crv.fit_iv(1);
			p0 = frms{in0}.End';
			p1 = frms{in1}.Start';
			leftout = in1 - in0 - 1;
			crv1 = [];
			crv1.ivid = 0;
			if leftout < 1 %% no gap
				crv1.iv = 0; crv1.fit_iv = 0;
				crv1.coeff = FIT.p2coeff(p0,p1,0,1);
				[crv1.Fit,~,crv1.len] = myTrajRender(size2(im), crv1.coeff, [0 1]);
				crv1.curve = p0;
				crv1.type='connect';
			else
				crv1.iv = [in0+1 in1-1];
				crv1.fit_iv = crv1.iv;
				pnt = FIT.find_intersect(crv0.coeff, crv.coeff, crv1.fit_iv, im);
				crv1.coeff = FIT.p3coeff(p0,pnt',p1);
				[crv1.Fit,~,crv1.len] = myTrajRender(size2(im), crv1.coeff, [0 1]);
				crv1.curve = pnt';
				crv1.type='bounce';
			end
			curves = [curves(1:end-1) crv1 curves(end)];
		end
		if kki > 2 %% connect between one serve
			crv0 = curves(end-1);
			in0 = crv0.fit_iv(2);
			in1 = crv.fit_iv(1);
			p0 = frms{in0}.End';
			p1 = frms{in1}.Start';
			p = curve(:,break_points(kki-1));
			leftout = in1 - in0 - 1;
			crv1 = [];
			crv1.iv = crv.iv;
			crv1.ivid = crv.ivid;
			remove_end = false;
			remove_end1 = false;
			if leftout == 1 %% only one unexplained
				in = crv0.fit_iv(2) + 1;
				crv1.fit_iv = [in in];
				crv1.coeff = FIT.p3coeff(p0,p,p1);
				[crv1.Fit,~,crv1.len] = myTrajRender(size2(im), crv1.coeff, [0 1]);
			elseif leftout == 0 %% no unexplained -> make bounce
				in0 = crv0.fit_iv(2);
				in1 = crv.fit_iv(1);
				p0 = frms{in0}.End';
				p1 = frms{in1}.Start';
				p = curve(:,break_points(kki-1));
				if norm(p - p0) < norm(p-p1) %% bounce at crv0
					curves(end-1).fit_iv(2) = curves(end-1).fit_iv(2)-1;
					if curves(end-1).fit_iv(1) > curves(end-1).fit_iv(2) 
						remove_end1 = true;
					else
						st_fit = curves(end-1).fit_iv(1); 
						en_fit = curves(end-1).fit_iv(2);
						max_power = max(2, min(6, round((en_fit-st_fit)/power_factor) ));
						if en_fit-st_fit == 0, max_power = 1; end
						curves(end-1).coeff = FIT.lsq_fit(frms, st_fit, en_fit, max_power, expos);
						[curves(end-1).Fit,~,curves(end-1).len] = myTrajRender(size2(im), curves(end-1).coeff, ...
								[curves(end-1).fit_iv(1) curves(end-1).fit_iv(2)+1]);
					end
					if in0 > 1 && ~isempty(frms{in0-1})
						p0 = frms{in0-1}.End';
					else
						p0 = frms{in0}.Start';
					end
					crv1.fit_iv = [curves(end-1).fit_iv(2)+1 crv.fit_iv(1)-1];
				else %% bounce at crv
					curves(end).fit_iv(1) = curves(end).fit_iv(1)+1;
					if curves(end).fit_iv(1) > curves(end).fit_iv(2) 
						remove_end = true;
					else
						st_fit = curves(end).fit_iv(1); 
						en_fit = curves(end).fit_iv(2);
						max_power = max(2, min(6, round((en_fit-st_fit)/power_factor) ));
						if en_fit-st_fit == 0, max_power = 1; end
						curves(end).coeff = FIT.lsq_fit(frms, st_fit, en_fit, max_power, expos);
						[curves(end).Fit,~,curves(end).len] = myTrajRender(size2(im), curves(end).coeff, ...
								[curves(end).fit_iv(1) curves(end).fit_iv(2)+1]);
					end
					if in1+1 <= numel(frms) && ~isempty(frms{in1+1})
						p1 = frms{in1+1}.Start';
					else
						p1 = frms{in1}.End';
					end
					crv1.fit_iv = [curves(end-1).fit_iv(2)+1 curves(end).fit_iv(1)-1];
				end
				crv1.coeff = FIT.p3coeff(p0,p,p1);
				[crv1.Fit,~,crv1.len] = myTrajRender(size2(im), crv1.coeff, [0 1]);
			else
				crv1.fit_iv = [in0+1 in1-1];
				p = FIT.find_intersect(crv0.coeff, crv.coeff, crv1.fit_iv, im)';
				crv1.coeff = FIT.p3coeff(p0,p,p1);
				[crv1.Fit,~,crv1.len] = myTrajRender(size2(im), crv1.coeff, [0 1]);
			end
			crv1.curve = p;
			crv1.type='bounce';
			crv_before = curves(1:end-1);
			crv_after = curves(end);
			if remove_end
				crv_after = [];
			end
			if remove_end1
				crv_before = curves(1:end-2);
			end
			curves = [crv_before crv1 crv_after];
		end
		%% predict to future
		if kki == numel(break_points) && ivid == numel(intervals)
			if curves(end).fit_iv(2) < numel(frms)
				crv = [];
				crv.iv = curves(end).iv;
				crv.ivid = curves(end).ivid;
				crv.coeff = curves(end).coeff;
				crv.fit_iv = [curves(end).fit_iv(2)+1 numel(frms)];
				if crv.fit_iv(2) - crv.fit_iv(1) > 20, continue; end
				crv.curve = curves(end).curve;
				[crv.Fit,~,crv.len] = myTrajRender(size2(im), crv.coeff, [crv.fit_iv(1) crv.fit_iv(2)+1]);
				crv.type = 'prediction';
				curves = [curves crv];
			end 
		end

	end
end



