function [x,y,z]=q_to_xyz(q,table)
    x=floor(table(q)/1000000);
    y=floor((table(q)-x*1000000)/1000);
    z=table(q)-x*1000000-y*1000;        
return;