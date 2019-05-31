function Data = generateData( nLine, noiseStd, rateOutlier,Uc)

	% parameters
    width = 640;
    height = 480;
	focal = 800;
	noiseStd   = noiseStd / focal; 

	% generate Data
% 	for kk= 1:nTest		    
		
		% lines
		[p1 P1_c] = randPts(nLine, width, height, focal, noiseStd,Uc);
		[p2 P2_c] = randPts(nLine, width, height, focal, noiseStd,Uc);
				
		%Rotation matrix from camera to world frame R_cw
		R_cw = randR();
		R_wc = R_cw.';
		
		% translation from camera to world frame T_cw
		T_cw = mean([P1_c P2_c],2);
		T_wc = -R_wc * T_cw;
				
		% world frame
		P1_w = R_wc * P1_c + kron(ones(1,nLine), T_wc);
		P2_w = R_wc * P2_c + kron(ones(1,nLine), T_wc);	

		% add outliers
		if rateOutlier > 0
			[i1,i2] = randOutlier(nLine,rateOutlier);
			p1(:,i1) = p1(:,i2);
			p2(:,i1) = p2(:,i2);
		end
		
		% save data
		Data.p1 = p1;
		Data.p2 = p2;
		Data.P1_w = P1_w;
		Data.P2_w = P2_w;
		Data.R_cw = R_cw;
		Data.T_cw = T_cw;
		
%     end
end

function [pt,PT] = randPts(n,w,h,f,noiseStd,Uc)


	if Uc==0 % uncentred
		u = rand(1,n) * w * 0.25;
		v = rand(1,n) * h * 0.25;
    elseif Uc==1 % centred
		u = rand(1,n) * w;
		v = rand(1,n) * h;
	end


	% depth range [4,8]
	depth = 4*rand(1,n)+4;
	
	% normalize
	u = (u - (w*0.5)) / f;
	v = (v - (h*0.5)) / f;

	% 2d
	pt = [u; v];

	% 3d
	PT = [u; v; ones(1,n)] .* kron(ones(3,1),depth);

	%add some noise to the endpoints
	if noiseStd > 0
		pt = pt + noiseStd*(rand(2,n)-0.5);
    end
end
		
function [i1,i2] = randOutlier(n,rateOutlier)

	m = round(n * rateOutlier);
	i0 = randperm(n);
	i1 = i0(1:m);
	i2 = i1(randperm(m));
end

function R = randR

    theta_x = 2*rand(1)-1;
    theta_y = 2*rand(1)-1;
    theta_z = 2*rand(1)-1;
    r_x = [1,0,0; 0, cos(theta_x), - sin(theta_x); 0, sin(theta_x), cos(theta_x) ];
    r_y = [cos(theta_y), 0, sin(theta_y); 0, 1, 0; - sin(theta_y), 0, cos(theta_y)];
    r_z = [cos(theta_z), - sin(theta_z), 0; sin(theta_z), cos(theta_z), 0; 0, 0, 1];
    R = r_z * r_y * r_x;
end