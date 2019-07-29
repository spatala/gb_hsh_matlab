function PI_ab = get_PIab(a,b)
% if a == b
%     PI_ab = sqrt(2*a+1);
% else
%     if (a > b)
%         PI_ab = sqrt(prod(2*b+1:2*a+1));
%     else
%         PI_ab = sqrt(prod(2*a+1:2*b+1));
%     end
% end

PI_ab = (sqrt(2*a+1))*(sqrt(2*b+1));

end