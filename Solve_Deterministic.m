function Solve_Deterministic(process)
    [x,y,theta,index,point,center,R] = Import_Data(process);
    G = zeros(length(x),1);
    for i = 1:length(x)
        supp_point(1) = x(i); supp_point(2) = y(i);
        G(i) = Pseudospectral_G(point,length(x), theta(i), center, R);
        disp(i)
    end
    B = Pseudospectral_B(point, center, R);
    file_output_G = append(append('Matlab_buffer/G_',string(process)),'.txt');
    fileID = fopen(file_output_G,'w');
    fprintf(fileID,'%f\n',G);
    fclose(fileID);
    file_output_B = append(append('Matlab_buffer/B_',string(process)),'.txt');
    fileID = fopen(file_output_B,'w');
    fprintf(fileID,'%f\n',B);
    fclose(fileID);
return