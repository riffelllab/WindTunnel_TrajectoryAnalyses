

%Plot the total amount of insects detected in each plot (as a bar plot)
function insects_in_frame(attr_frame)
    figure()
    nbins=[min(attr_frame):max(attr_frame)];
    hist(double(attr_frame), double(nbins))
    title('Number of insect detected in each frame')
    xlabel('Frame number');
    ylabel('total number of insect detected');
end

y = zeros(size(attr_frame));
for i = 1:length(attr_frame)
y(i) = sum(attr_frame==attr_frame(i));
end
plot(y)