A = load("result-dense-discrete.txt");
x1 = cos(A(:,2));
y1 = -sin(A(:,2));
x2 = x1 + 0.8*cos(A(:,3));
y2 = y1 - 0.8*sin(A(:,3));
figure();
[A,map] = rgb2ind(frame2im(getframe),256);
imwrite(A,map,'1.gif','LoopCount',65535,'DelayTime',0.1);
for ii = 1:55
    clf;
    h = plot(x2(1:ii), y2(1:ii));
    hold on
    axis equal
    axis([-2,2,-2,2]);
    l1 = line([0,x1(ii)], [0,y1(ii)]);
    l2 = line([x1(ii),x2(ii)], [y1(ii),y2(ii)]);
    [A,map] = rgb2ind(frame2im(getframe),256);
    imwrite(A,map,'1.gif','WriteMode','append','DelayTime',0.02);
end
for ii = 56:length(x2)
    clf;
    h = plot(x2(ii-55:ii), y2(ii-55:ii));
    hold on
    axis equal
    axis([-2,2,-2,2]);
    l1 = line([0,x1(ii)], [0,y1(ii)]);
    l2 = line([x1(ii),x2(ii)], [y1(ii),y2(ii)]);
    [A,map] = rgb2ind(frame2im(getframe),256);
    imwrite(A,map,'1.gif','WriteMode','append','DelayTime',0.02);
end
