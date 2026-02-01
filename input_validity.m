function [valid]=input_validity( ap, Delta_x, s )
%imposes restriction on the input
    valid=true;
    for n=1:s+1
        if ap(1,n)<Delta_x
            valid=false;
            return
        end
    end
    if ap(1,s+1) > 1 || ap(1,s+1) < 0.1
        valid=false;
        return
    end
end
