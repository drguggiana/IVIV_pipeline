function [alpha_ang]=sign_ang(data_in)

for m=1:length(data_in)
              if data_in(m,10)==1 
                  temp(:,m)=(90-abs(data_in(m,5)))-180
              elseif data_in(m,10)==2 
                  temp(:,m)=abs((90-abs(data_in(m,5)))-180);
                   elseif data_in(m,10)==3 
                       temp(:,m)=(90-abs(data_in(m,5)))*-1
                        elseif data_in(m,10)==4 
                            temp(:,m)=(90-abs(data_in(m,5)))
              end
end
   
alpha_ang=temp';
end