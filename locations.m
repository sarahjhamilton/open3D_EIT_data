points=[4.25	-12.75	4.25
-4.25	-12.75	4.25
-4.25	-12.75	-4.25
4.25	-12.75	-4.25
4.25	-8.5	8.5
-4.25	-8.5	8.5
-8.5	-8.5	4.25
-8.5	-8.5	-4.25
-4.25	-8.5	-8.5
4.25	-8.5	-8.5
8.5	-8.5	-4.5
8.5	-8.5	4.5
4.25	0	8.5
-4.25	0	8.5
-8.5	0	4.25
-8.5	0	-4.25
-4.25	0	-8.5
4.25	0	-8.5
8.5	0	-4.5
8.5	0	4.5
4.25	8.5	8.5
-4.25	8.5	8.5
-8.5	8.5	4.25
-8.5	8.5	-4.25
-4.25	8.5	-8.5
4.25	8.5	-8.5
8.5	8.5	-4.5
8.5	8.5	4.5
4.25	12.75	4.25
-4.25	12.75	4.25
-4.25	12.75	-4.25
4.25	12.75	-4.25];

leftend=[-8.5 -12.75 -8.5
    -8.5 -12.75 8.5
    8.5 -12.75 8.5
    8.5 -12.75 -8.5
    -8.5 -12.75 -8.5];

rightend=[-8.5 12.75 -8.5
    -8.5 12.75 8.5
    8.5 12.75 8.5
    8.5 12.75 -8.5
    -8.5 12.75 -8.5];

side1=[-8.5 -12.75 8.5
    -8.5 12.75 8.5];
side2=[-8.5 -12.75 -8.5
    -8.5 12.75 -8.5];
side3=[8.5 -12.75 8.5
    8.5 12.75 8.5];
side4=[8.5 -12.75 -8.5
    8.5 12.75 -8.5];

figure;
for i=1:4
    plot3(points(i,1),points(i,2),points(i,3),'or')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=29:32
    plot3(points(i,1),points(i,2),points(i,3),'or')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=5:8:21
    plot3(points(i,1),points(i,2),points(i,3),'*c')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=6:8:22
    plot3(points(i,1),points(i,2),points(i,3),'*c')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=7:8:23
    plot3(points(i,1),points(i,2),points(i,3),'^k')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end

for i=8:8:24
    plot3(points(i,1),points(i,2),points(i,3),'^k')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=9:8:29
    plot3(points(i,1),points(i,2),points(i,3),'sm')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=10:8:26
    plot3(points(i,1),points(i,2),points(i,3),'sm')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=11:8:27
    plot3(points(i,1),points(i,2),points(i,3),'pb')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end
for i=12:8:29
    plot3(points(i,1),points(i,2),points(i,3),'pb')
    txt=strcat({'  '},num2str(i));
    text(points(i,1),points(i,2),points(i,3),txt)
    hold on
end


plot3(leftend(:,1),leftend(:,2),leftend(:,3))
plot3(rightend(:,1),rightend(:,2),rightend(:,3))

plot3(side1(:,1),side1(:,2),side1(:,3))
plot3(side2(:,1),side2(:,2),side2(:,3))
plot3(side3(:,1),side3(:,2),side3(:,3))
plot3(side4(:,1),side4(:,2),side4(:,3));
axis('equal')
