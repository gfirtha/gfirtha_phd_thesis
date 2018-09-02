function [ dref,ref_win_x0 ] = find_ref_dist( x_ssd, reference_curve, K )
wb = waitbar(0,'Calculating');

TA = zeros(2,size(reference_curve,1)-1);
dref = zeros(size(x_ssd,1),1);
ref_win_x0  = zeros(size(x_ssd,1),1);
for n = 1 : length(x_ssd)
    waitbar(n/length(x_ssd),wb);
    
    x0 = x_ssd(n,:)';
    v =  K(n,:)';
    
    for m = 1 : size(reference_curve,1)-1
        x1 = reference_curve(m,:)';
        x2 = reference_curve(m+1,:)';
        A = [ v x2-x1 ];
        TA(:,m) = A\(x2-x0);
    end
    d = TA(1,:);
    a = TA(2,:);
    d( a<0 | a>1 ) = [];
    a( a<0 | a>1 ) = [];
    if ~isempty(d)
        dref(n) = min(d);
        ref_win_x0(n)  = 1;
    end
end
close(wb);

end

