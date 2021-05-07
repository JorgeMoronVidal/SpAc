function [x,y,theta,index,point,center,R] = Import_Data(process)
    delimiterIn = ',';
    headerlinesIn = 0;
    parameters_route = append(append('Matlab_buffer/parameters_',string(process)),'.csv');
    parameters_table = importdata(parameters_route,delimiterIn, headerlinesIn);
    point(1) = parameters_table(01,1);
    point(2) = parameters_table(1,2);
    center(1) = parameters_table(1,3);
    center(2) = parameters_table(1,4);
    R = parameters_table(1,5);
    points_route = append(append('Matlab_buffer/points_',string(process)),'.csv');
    points_table = importdata(points_route,delimiterIn, headerlinesIn);
    index = points_table(:,1);
    theta = points_table(:,2);
    x = center(1) + cos(theta);
    y = center(2) + sin(theta);
return 