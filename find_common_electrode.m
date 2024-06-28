function list = find_common_electrode(position1, position2)
N = size(position1,1);

count = 0;
list = zeros(1,2);
for i = 1:N
    diff_pos = sum(abs(position2 - position1(i,:)),2);
    j = find(diff_pos<eps);
    if ~isempty(j)
        count = count + 1;
        list(count,1) = i;
        list(count,2) = j;
    end
end