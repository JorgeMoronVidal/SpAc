function [G_i,G_j,G_ij] = Print_Dependencies(Circle_position, N_corner, N_side, N_interior, N_circles, x,y,G,G_i,G_j,G_ij)
    [i_begin,i_end] = Compute_knot_indexes(Circle_position, N_corner, N_side, N_interior, N_circles);
    G_j_patt = [i_begin:i_end];
    figure
    hold on
    scatter(x(i_begin:i_end),y(i_begin:i_end),'filled','b')
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1)-1,Circle_position(2)-1], N_corner, N_side, N_interior, N_circles);
    scatter(x(i_begin+1:i_begin+8),y(i_begin+1:i_begin+8),'filled','r')
    for i = i_begin+1:i_begin+8
        G_i = [G_i; i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1)-1,Circle_position(2)], N_corner, N_side, N_interior, N_circles);
    scatter(x(i_begin:i_begin+5),y(i_begin:i_begin+5),'filled','r')
    for i = i_begin:i_begin+5
        G_i = [G_i;i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    scatter(x(i_end-6:i_end),y(i_end-6:i_end),'filled','r')
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1)-1,Circle_position(2)+1], N_corner, N_side, N_interior, N_circles);
    for i = i_end-6:i_end
        G_i = [G_i; i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    scatter(x(i_end-9:i_end-2),y(i_end-9:i_end-2),'filled','r')
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1),Circle_position(2)+1], N_corner, N_side, N_interior, N_circles);
    for i = i_begin + 26:i_begin+38
        G_i = [G_i; i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    scatter(x(i_begin + 26:i_begin+38),y(i_begin + 26:i_begin+38),'filled','r')
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1)+1,Circle_position(2)+1], N_corner, N_side, N_interior, N_circles);
    for i = i_end-20:i_end-13
        G_i = [G_i; i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    scatter(x(i_end-20:i_end-13),y(i_end-20:i_end-13),'filled','r')
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1)+1,Circle_position(2)], N_corner, N_side, N_interior, N_circles);
    for i = i_begin+15:i_begin+27
        G_i = [G_i; i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    scatter(x(i_begin+15:i_begin+27),y(i_begin+15:i_begin+27),'filled','r')
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1)+1,Circle_position(2)-1], N_corner, N_side, N_interior, N_circles);
    for i = i_begin+12:i_begin+19
        G_i = [G_i; i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    scatter(x(i_begin+12:i_begin+19),y(i_begin+12:i_begin+19),'filled','r')
    [i_begin,i_end] = Compute_knot_indexes([Circle_position(1),Circle_position(2)-1], N_corner, N_side, N_interior, N_circles);
    for i = i_begin+4:i_begin+16
        G_i = [G_i; i];
        G_j = [G_j;i];
        G_ij = [G_ij;1.0];
        for j = G_j_patt
            G_i = [G_i; i];
            G_j = [G_j;j];
            G_ij = [G_ij;1.0];
        end
    end
    scatter(x(i_begin+4:i_begin+16),y(i_begin+4:i_begin+16),'filled','r')
    figure
    G_aux = sparse(G_i,G_j,G_ij);
    spy(G_aux);
end