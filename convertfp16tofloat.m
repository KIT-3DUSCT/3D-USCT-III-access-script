function output=convertfp16tofloat(input)

output=zeros(size(input));
i2=input;
input=typecast(input(:),'uint16');
%%%%%%%%%%%%%%% Looped SLOW
% for j=1:size(input,2)
%     %ascan=input(:,j);
%     j
%     for i=1:size(input,1)
%         fp16=input(i,j);
%       %fp16=typecast(fp16,'uint16');
% 
%         if fp16>=2^15 sign_bit=1; else sign_bit=0; end
%         exponent = bitand(bitshift(fp16,-11),15); %  4 bits
%         fraction = bitand(fp16 , 2047); % 11 bits
%         hidden_bit11 = uint16(0);
% 
%         if (exponent == 0) hidden_bit11 = uint16(0); % value too small
%            if (sign_bit == 0)
%                 hidden_bit11 = uint16(0); % positive
%             else
%                 hidden_bit11 = 2048; % negative
%             end
%         else
%             if (sign_bit == 0)
%                 hidden_bit11 =  uint16(2048); % positive
%             else hidden_bit11 = uint16(0); % negative
%             end
%         end
% 
%         fraction= fraction + hidden_bit11; % now we’ve got 12bits
%         
%         if (exponent == 0) % value too small
%            if (sign_bit == 0) % positive
%                 Z = double(fraction);
%             else % negative
%                 %Z =-double(fraction);
%                % Z = double((-1)*bitand( ((int32(2^16)-int32(fraction)) + 1) ,  int32(4095) ));
%                 Z = double( int32(fraction)-2^(11));
%             end
%         else % exponent > 0
%             if (sign_bit == 0) % positive
%                 Z = double(fraction) * 2^(double(exponent)-1);
%             else %% negative
%                 %Z =-double(fraction)* 2^(double(exponent)-1);
%                 %Z = double((-1) * bitand(uint32((2^16-int32(fraction)) + 1) ,  int32(4095)) * 2^int32(exponent-1));
%                 Z = double( (int32(fraction)-2^(11)) * 2^(int32(exponent)-1) );
%             end
%         end
% 
%         %%write back
%         output(i,j)=Z;
% 
%     end
% end %%%loop over data


%%%%%%%%%%%%%%%%%vectorized
  exponent = bitand(bitshift(typecast(input,'uint16'),repmat(-11,size(input))),repmat(uint16(15),size(input))); %  4 bits
  sign_bit = double(typecast(input,'int16')<0);
  %fraction = bitand(typecast(input,'uint16') , 2047); % 11 bits
  %fraction2 = fix(double(typecast(input,'int16'))./(2.^5)); % 11 bits
  fraction3 = uint32(bitand(typecast(input,'uint16'),repmat(uint16(2047),size(input))));
  hidden_bit = uint32((sign_bit==1).*(exponent==0).*2048 + (sign_bit==0 & exponent~=0).*2048); %(-sign_bit+(sign_bit==0)).*(-(exponent==0)+(exponent~=0)).*2048;
  output=typecast(bitshift((fraction3 + hidden_bit + uint32(repmat((2^32-1)-sum(2.^(0:11)),size(input)).*sign_bit)),uint32(exponent-1)),'int32'); % + uint32(2.^double((exponent+1)).*sign_bit) 
  %%%% 
  output=reshape(output,size(i2));
  %%%% fraction3 = double(mod(typecast(input,'int16'),int16((-sign_bit+(sign_bit==0))*2^10)));
  %%%% hidden_bit =(-sign_bit+(sign_bit==0)).*(exponent~=0).*2048 - (sign_bit==1 & exponent==0)*4096; %(-sign_bit+(sign_bit==0)).*(-(exponent==0)+(exponent~=0)).*2048;
  %%%% output=(fraction3+hidden_bit).*2.^double(exponent-1);
end

