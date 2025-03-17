function plotMesh( file )

h = figure(1);
hold on
for e = 1: file.no_of_elements
   for k = 1 : size(file.data{e},1)
       for j = 1 : size(file.data{e},2)
           for i = 1 : size(file.data{e},3)
               x = file.data{e}(k,j,i,:);
               plot3(x(1),x(2),x(3),'.k');
           end
       end
   end
end
    
end

