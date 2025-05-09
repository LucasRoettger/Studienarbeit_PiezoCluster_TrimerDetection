function M = getRotationMatrix(alfa, axes)
    switch axes
        case 'X'
            M = [1 0 0; 0 cos(alfa) -sin(alfa); 0 sin(alfa) cos(alfa)];
        case 'Y'
            M = [cos(alfa) 0 sin(alfa); 0 1 0; -sin(alfa) 0 cos(alfa)];
        case 'Z'
            M = [cos(alfa) -sin(alfa) 0; sin(alfa) cos(alfa) 0; 0 0 1];
        otherwise
            error("'axes' statement must be 'X', 'Y', or'Z'!")
    end
end