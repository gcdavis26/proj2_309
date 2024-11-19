function [dcm] = dcm3(angle)
dcm = [cosd(angle),sind(angle),0;-sind(angle),cosd(angle),0;0,0,1];
end