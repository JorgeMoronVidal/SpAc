function [knot_index_begin,knot_index_end] = Compute_knot_indexes(Circle_position, N_corner, N_side, N_interior, N_circles)
    if Circle_position(1) == 0
        if(Circle_position(2) == 0)
            knot_index_begin = 1;
            knot_index_end = N_corner;
        else
            if(Circle_position(2) == N_circles(2) -1)
                knot_index_begin = N_corner + (N_circles(2)-2)*N_side + 1;
                knot_index_end = 2*N_corner + (N_circles(2)-2)*N_side;
            else
                knot_index_begin = N_corner + (Circle_position(2)-1)*N_side + 1;
                knot_index_end = N_corner + (Circle_position(2))*N_side;
            end
        end
    else
        if Circle_position(1) == N_circles(1) -1
            if(Circle_position(2) == 0)
                knot_index_begin = 2*N_corner + (N_circles(2)-2 + 2*(N_circles(1)-2))*N_side + (N_circles(2)-2)*(N_circles(1)-2)*N_interior +1;
                knot_index_end = 3*N_corner + (N_circles(2)-2 + 2*(N_circles(1)-2))*N_side + (N_circles(2)-2)*(N_circles(1)-2)*N_interior;
            else
                if(Circle_position(2) == N_circles(2) -1)
                    knot_index_begin = 3*N_corner + (2*(N_circles(2)-2) + 2*(N_circles(1)-2))*N_side + (N_circles(2)-2)*(N_circles(1)-2)*N_interior + 1;
                    knot_index_end = 4*N_corner + (2*(N_circles(2)-2) + 2*(N_circles(1)-2))*N_side + (N_circles(2)-2)*(N_circles(1)-2)*N_interior;
                else
                    knot_index_begin = 3*N_corner + (2*(N_circles(2)-2) + (N_circles(1)-2) + Circle_position(2)-1)*N_side + (N_circles(2)-2)*(N_circles(1)-2)*N_interior + 1;
                    knot_index_end = 3*N_corner + (2*(N_circles(2)-2) + (N_circles(1)-2) + Circle_position(2))*N_side + (N_circles(2)-2)*(N_circles(1)-2)*N_interior;
                end
            end        
        else
            if(Circle_position(2) == 0)
                knot_index_begin = 2*N_corner + (N_circles(2)-2 + 2*(Circle_position(1)-1))*N_side +(N_circles(2)-2)*(Circle_position(1)-1)*N_interior +1;
                knot_index_end = 2*N_corner + (N_circles(2)-2 + 2*(Circle_position(1)-1) + 1)*N_side +(N_circles(2)-2)*(Circle_position(1)-1)*N_interior;
            else
                if(Circle_position(2) == N_circles(2) -1)
                    knot_index_begin = 2*N_corner + (N_circles(2)-2 + 2*(Circle_position(1)-1) + 1)*N_side +(N_circles(2)-2)*(Circle_position(1))*N_interior +1;
                    knot_index_end = 2*N_corner + (N_circles(2)-2 + 2*(Circle_position(1)))*N_side +(N_circles(2)-2)*(Circle_position(1))*N_interior;
                else
                    knot_index_begin = 2*N_corner + (N_circles(2)-2 + 2*(Circle_position(1)-1) + 1)*N_side +((N_circles(2)-2)*(Circle_position(1)-1)+(Circle_position(2)-1))*N_interior +1;
                    knot_index_end = 2*N_corner + (N_circles(2)-2 + 2*(Circle_position(1)-1) + 1)*N_side +((N_circles(2)-2)*(Circle_position(1)-1)+(Circle_position(2)))*N_interior;
                end
            end
        end
    end
end