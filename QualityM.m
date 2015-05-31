

% --------------------------------------------------------------------------------------------------
%
%    this code is a part of the demo software for Performance comparison of the ACFs produced by the LFM chirp windowed by the 
%    Taylor window, and the optimized NLFM signal ,to reproduce some results described in the papers:
%
%    Zhi-huo XU, Yun-kai Deng, Robert Wang, et al., "A New Optimization Framework for Nonlinear FM Chirps Using the Firefly Algorithm", 
%
%    submitted to IEEE Geosci. Remote Sens. Lett.
%
% --------------------------------------------------------------------------------------------
%
% authors:              Zhi-huo Xu, Yun-kai Deng,Robert Wang, et al.
%
% web page:              http://
%
% contact:               xuzhihuo@gmail.com
%
% --------------------------------------------------------------------------------------------
% Copyright (c) 2015 IECAS.
% Institute of Electronics, Chinese Academy of Sciences.
% All rights reserved.
% This work should be used for nonprofit purposes only.
% --------------------------------------------------------------------------------------------


function  [pslr, idx_l,idx_r]=QualityM(amp)

n=length(amp);

idx(n)=0;
for k=2:n-1
    
    v_mid=amp(k);
    v_left=amp(k-1);
    v_right=amp(k+1);
    
    if(v_mid<v_left)&&(v_mid<v_right)
        idx(k)=1;
    end
end


idx_mid=find(amp==max(amp));
idx_mid=idx_mid(1);   

for k=idx_mid:n-1
    
    if(1==idx(k))
        idx_r=k;
        break;
    end
end

for k=idx_mid:-1:1
    if(1==idx(k))
        idx_l=k;
        break;
    end
end

amp(idx_l:idx_r)=-147;

idx_mid=find(amp==max(amp));

pslr=amp(idx_mid(1));
