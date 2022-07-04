function m=makehgtform(varargin)
%makehgtform('zrotate',2*pi*(rotshift+movement(1))/360)*makehgtform('translate',[0 0 movement(2)]);
   m = eye(4);
   
if strcmp(varargin{1},'translate')
                y = varargin{2}(2);
                z = varargin{2}(3);
                x = varargin{2}(1);
                    m = m * [1 0 0 x; ...
                     0 1 0 y; ...
                     0 0 1 z; ...
                     0 0 0 1];
                    
 end
 
 if strcmp(varargin{1},'zrotate')
  t = varargin{2};
     ct = cos(t);
            st = sin(t);
            m = m * [ct -st 0 0; ...
                     st  ct 0 0; ...
                      0   0 1 0; ...
                      0   0 0 1];
end
  
end