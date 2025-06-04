close all; clear all;

plot = figure("inverthardcopy", "off");

low_index = 1;
high_index = 50;

  file_path = 'output/2_output.dat';
  header_lines = 1;  % Number of header lines to skip
  % Read the data skipping the header line
  output = dlmread(file_path, '', header_lines, 0);

  beta = sqrt(output(:,6).^2+output(:,7).^2+output(:,8).^2)/2.997925e8;
  color_vals = ones(size(output(:,1))) .* beta;

  if (length(output(:,1)) > 1)
    [sx sy sz] = sphere(20);
    sphere_colors = ones(size(sz))*beta(1);
    surf([output(:,3) output(:,3)], [output(:,4) output(:,4)], [output(:,5) output(:,5)], [color_vals color_vals], "FaceColor", "none", "EdgeColor", "interp", "linewidth", 2);
    hold on;
    surf(sx*0.25 + output(1,3), sy*0.25 + output(1,4), sz*0.25 + output(1,5), sphere_colors, "EdgeColor", "none");
  endif
  hold on;

 L_id = fopen("../Magnetic_Code/Laplacian.in", "r");
 Laplacian_data = textscan(L_id, "%f", "CommentStyle", "!");
 Laplacian_data = cell2mat(Laplacian_data); % convert cell array to numeric array
 fclose(L_id);
 I_id = fopen("input.dat","r");
 input_data = textscan(I_id, "%f", "CommentStyle", "//"); %reads in everything as floats
 input_data = cell2mat(input_data); % convert cell array to numeric array
 fclose(I_id);;

 grid_size = Laplacian_data(1);
 factor = Laplacian_data(2);
 grid_SI = grid_size*factor;
 cyl_rad = factor*Laplacian_data(4);
 cyl_hei = 2*factor*Laplacian_data(7);

 input_data(1) = round(input_data(1)); % convert number of runs to an int

 %compute where the run info starts
 input_line_factor = input_data(1)*6+1 % number of particle values plus first row

 sph_rad = input_data(input_line_factor+6);
 disp(sph_rad);


 center = grid_SI/2;

 [cx cy cz] = cylinder([cyl_rad cyl_rad], 50);
 [sx sy sz] = sphere(20);
 circx = cx; circx(1,:) = 0;
 circy = cy; circy(1,:) = 0;
 circzt = cz * 0;
 circzb = cz * 0;

% legend("show", "location", "eastoutside");

 set(gca,"fontsize", 16, "color", [0.77, 0.77, 0.77]);

 %COLORMAP
 % Define the three colors for the gradient
color1 = [1, 0, 1];       % ?
color2 = [0, 0, 1];       % Blue
color3 = [1, 0, 0];       % Red

% Number of colors in the colormap
num_colors = 256;
% Create a custom colormap with the gradient
custom_colormap = interp1([1, num_colors/2, num_colors], [color1; color2; color3], 1:num_colors);
 colormap("hot");
% colormap("hot");
 caxis([0 1]);
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
 view(-30, 15);
 hold off;
 grid on;
 print(plot, "output_singleRun.png", "-dpng");
 pause;
