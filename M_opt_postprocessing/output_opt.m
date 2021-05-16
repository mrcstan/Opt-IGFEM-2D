%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Output2 Module%%%%%%%%%%%%%%%%%%%%%%

function [] = output_opt (param_val, theta, gradTheta, opt_outfile,counter)

fileID = fopen(opt_outfile, 'w');

%%%%%%% part (e)-i
fprintf (fileID, '\n\n*************************************************');
fprintf (fileID, '\n*************** Optimization ********************');
fprintf (fileID, '\n*************************************************');
% fprintf (fileID, ' \n(d) Optimize the beam shape to minimize the compliance subject to ');
% fprintf (fileID, ' the constraint that the area should not exceed\n ');
% fprintf (fileID, ' l × w = 2m2 . Denote the compliance of the initial and optimized designs.\n\n');
fprintf(fileID, '\n');

fprintf(fileID, ' Function evaluation: %i \n\n',  counter);

fprintf (fileID, ' param_val\n ');
for i=1:numel(param_val)
    fprintf(fileID, ' %16.15e \n', param_val(i));
end
fprintf(fileID, '\n');

% % param_name = {'a: ', 'b: ', 'c: ', 'Xc: ', 'Yc: '};
% for i=1:size(interface, 2)
%     fprintf (fileID, ' a: ');
%     fprintf(fileID, ' %16.15e \n', ...
%                   interface(1, i).constants.a);
%     fprintf (fileID, ' b: ');
%     fprintf(fileID, ' %16.15e \n', ...
%                   interface(1, i).constants.b);
%     fprintf (fileID, ' c: ');
%     fprintf(fileID, ' %16.15e \n', ...
%                   interface(1, i).constants.c);
%     fprintf (fileID, ' Xc: ');
%     fprintf(fileID, ' %16.15e \n', ...
%                   interface(1, i).constants.Xc);
%     fprintf (fileID, ' Yc: ');
%     fprintf(fileID, ' %16.15e \n', ...
%                   interface(1, i).constants.Yc);
% end
% fprintf(fileID, '\n');


fprintf (fileID, ' theta\n ');
for i=1:numel(theta)
    fprintf(fileID, ' %16.15e \n', theta(i));
end
fprintf(fileID, '\n');

fprintf (fileID, ' grad theta\n ');
for i=1:numel(gradTheta)
    fprintf(fileID, ' %16.15e \n', gradTheta(i));
end
fprintf(fileID, '\n');

% fprintf(fileID, ' grad_theta: \n');
% for i =1:n_funct
%     for j = 1:n_param
%         fprintf(fileID, ' %16.15e \t',grad_theta(i,j));
%     end
%     fprintf(fileID, '\n');
% end

fprintf(fileID, '\n');
fprintf (fileID, '\n*************************************************\n');
fprintf(fileID, '\n');

fclose (fileID);

return