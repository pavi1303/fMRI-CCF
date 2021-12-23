function [ Qdata,mask,table ] = xyz_to_q_ALT( data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    maxData = max(max(max(max(data))));
    [tdim,xres,yres,zres]=size(data);
    Qdata = zeros(tdim,xres*yres*zres);
    mask=zeros(xres,yres,zres);
    table=zeros(xres*yres*zres,1);
    q=0;
    for x=1:xres
        for y=1:yres
            for z=1:zres
                data_xyz = data(:,x,y,z);
                if mean(data_xyz)>0 && max(data_xyz)>0.1*maxData
                    q=q+1;
                    Qdata(:,q)=data(:,x,y,z);
                    mask(x,y,z)=q;
                    table(q,1)=x*1000000+y*1000+z;
                end
            end
        end
    end
    Qdata = Qdata(:,1:q);
end


  


