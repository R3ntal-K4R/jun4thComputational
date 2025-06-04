close all; clear all;

plot = figure("inverthardcopy", "off");

low_index = 1;
high_index = 229;

directory_range = low_index:high_index;
colors = [[0,0,0];[1,0,0];[0,0,1];[0,0.85,0];[136/255,46/255,114/255]; ...
  [247/255,203/255,69/255];[238/255,85/255,24/255];[114/155,25/155,14/255]; ...
  [119/255,119/255,119/255];[110/255,166/255,205/255]];

max_beta = 0;

for i = directory_range
  file_path = ["output/",num2str(i),'_output.dat'];
  disp(file_path);
%  disp(i);
%  disp(mod(i,11));
  header_lines = 1;  % Number of header lines to skip
  % Read the data skipping the header line
  output = dlmread(file_path, '', header_lines, 0);

  beta = sqrt(output(1,6).^2 + output(1,7).^2 + output(1,8).^2) ./ 2.997925e8;
  color_vals = ones(size(output(:,1))) * beta;

  % change max beta if applicable
  if (beta > max_beta)
    max_beta = beta;
  endif

  if (length(output(:,1)) > 1)
    [sx sy sz] = sphere(20);
    sphere_colors = ones(size(sz))*beta;
    surf([output(:,3) output(:,3)], [output(:,4) output(:,4)], [output(:,5) output(:,5)], [color_vals color_vals], "FaceColor", "none", "EdgeColor", "interp", "linewidth", 2);
    surf(sx*0.25 + output(1,3), sy*0.25 + output(1,4), sz*0.25 + output(1,5), sphere_colors, "EdgeColor", "none");
  endif
  hold on;

%  plot3(output(:,3),output(:,4),output(:,5),'-', "linewidth", 2, 'color', colors(mod(i,10)+1,:), "displayname", ['{\beta}=' num2str(beta)]);
%  hold on;
%  plot3(output(1,3),output(1,4),output(1,5), '*', "markersize", 8, 'color', colors(mod(i,10)+1,:));
 endfor

 L_id = fopen("../Magnetic_Code/Laplacian.in", "r");
 Laplacian_data = textscan(L_id, "%f", "CommentStyle", "!");
 Laplacian_data = cell2mat(Laplacian_data); % convert cell array to numeric array
 fclose(L_id);
 I_id = fopen("input.dat","r");
 input_data = textscan(I_id, "%f", "CommentStyle", "//"); %reads in everything as floats
 input_data = cell2mat(input_data); % convert cell array to numeric array
 fclose(I_id);;

 input_data(1) = round(input_data(1)); % convert number of runs to an int

 %compute where the run info starts
 input_line_factor = input_data(1)*6+1 % number of particle values plus first row

 grid_size = Laplacian_data(1);
 factor = Laplacian_data(2);
 grid_SI = grid_size*factor;
 cyl_rad = factor*Laplacian_data(4);
 cyl_hei = 2*factor*Laplacian_data(7);
 center = grid_SI/2;


 sph_rad = input_data(input_line_factor+6);
 % disp(sph_rad);

 [cx cy cz] = cylinder([cyl_rad cyl_rad], 50);
 [sx sy sz] = sphere(20);
 circx = cx; circx(1,:) = 0;
 circy = cy; circy(1,:) = 0;
 circzt = cz * 0;
 circzb = cz * 0;

% legend("show", "location", "eastoutside");

 set(gca,"fontsize", 16, "color", [0.77, 0.77, 0.77]);
 colormap("hot");
% caxis([0 1])        % color range as the entire beta range
 caxis([0 max_beta]); % have the color range be based on the max beta in the plotted files
 hcb = colorbar("location", "eastoutside", "fontsize", 16);
 colorTitleHandle = get(hcb,'Title');
 titleString = '{\beta}';
 set(colorTitleHandle ,'String',titleString, "fontsize", 16);
 cylinder_colors = ones(size(cz))*0.07;
 sphere_colors = ones(size(sz))*0.15;
 cap_colors = ones(size(circzt))*0.07;

 surf(cx+center, cy+center, cyl_hei*cz+center-cyl_hei/2, cylinder_colors); % Cylinder
 surf(circx+center, circy+center, circzt+center+cyl_hei/2, cap_colors); % Cylinder top
 surf(circx+center, circy+center, circzb+center-cyl_hei/2, cap_colors); % Cylinder bottom
 surf(sph_rad*sx+center, sph_rad*sy+center, sph_rad*sz+center+cyl_hei/2+sph_rad, sphere_colors, "FaceColor", "none"); % Top  sphere
 surf(sph_rad*sx+center, sph_rad*sy+center, sph_rad*sz+center-cyl_hei/2-sph_rad, sphere_colors, "FaceColor", "none"); % Bottom sphere

 axis([-grid_SI/4 grid_SI+grid_SI/4 -grid_SI/4 grid_SI+grid_SI/4 -grid_SI/4 grid_SI+grid_SI/4]);
 title("Particle Trajectories", "fontsize", 16);
% set(gca,"fontsize", 16, "color", [0.77, 0.77, 0.77]);
% view(-30, 15);  % diagonal view
 view(0, 0);      % side view
% view(0, 90);    % top view
 hold off;

 print(plot, "output_multiRuns.png", "-dpng");
 pause;
