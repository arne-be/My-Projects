%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              LABORATORY #6 
%%%              COMPUTER VISION 2021-2022
%%%              NON-RIGID STRUCTURE FROM MOTION - OPTIMIZATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [J]=JacobianPattern(K,n_frames,n_points,vij,priors)
% Input
% K: shape basis rank
% n_frames: number of frames 
% n_points: number of points 
% vij: visibility map
% - priors: structure with fields:
%         priors.camera_prior: boolean, 1 for rotation smoothness on, 0 off
%         priors.coeff_prior: boolean, 1 for deformation smoothness on, 0 off 
% Output
% J: the Jacobian matrix pattern

 
% J data term is shorter than 2xFxP if the FxP visibility contains zeros
prior_terms = priors.coeff_prior + priors.camera_prior;



% Prior_terms must be a number from 0 to 2
if prior_terms < 0 || prior_terms > 2
    error('wrong values in prior options');
end

% Jacobian matrix pattern definition
% I give you the size, but you need to define the "ones"
J = zeros(2*nnz(vij)+ prior_terms*(n_frames-1),(K+6)*n_frames + K*3*n_points);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%% MISSING CODE HERE %%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % You need to include the correct relations in J. First of all, consider the data term, and then, the priors
    
     % Data term
     % 

     %as info the output of model2vector
        % v: a vector with all parameters you should estimate, following the order
        % [{L_11, ..., L_K1, R_1 T_1} ... {L_1F, ..., L_KF, R_F T_F}, X_1, ..., X_K]
        % i.e., for every frame until F, to include weight coefficients rotations
        % (4 entries, quaternions) and translations (2 entries), and finally, the K shape vectors
    %as info the output of vector2model
        % v: a vector with all parameters you should estimate, following the order
        % [{L_11, ..., L_K1, R_1 T_1} ... {L_1F, ..., L_KF, R_F T_F}, X_1, ..., X_K]
        % i.e., for every frame until F, to include weight coefficients rotations
        % and translations, and finally, the K shape vectors


 for f=0:(n_frames-1)
    i_begin = f*(K+6)+1;
    j_begin = f*(2*n_points)+1;

    J(j_begin:j_begin + (2*n_points)-1, i_begin:i_begin + (K+6)-1) = 1;

    for k=0:(K-1)
        for i=0:2
            i_begin = (K+6)*n_frames + k*3*n_points + i*n_points + 1;
            j_begin = f*(2*n_points) + i*cast(n_points/3, "uint32") + 1;

            J(j_begin:j_begin + cast(n_points/3, "uint32") -1, i_begin:i_begin + n_points -1) = 1;
            
            J(j_begin+n_points:j_begin + n_points+cast(n_points/3, "uint32") -1, i_begin:i_begin + n_points -1) = 1;
        end
       
    end
 end
    
 if (priors.coeff_prior == 1)
     % prior terms on L
    coeff_p = zeros((n_frames-1)*2, (K+6)*n_frames + K*3*n_points);
    for f=0:n_frames-1
        i_begin = f*(K+6)+1;
        j_begin = f*K+1;

        coeff_p(j_begin:j_begin+1, i_begin:i_begin+1) = 1;
        coeff_p(j_begin:j_begin+1, i_begin+(6+K):i_begin+(6+K)+1) = 1;

    end
    J = cat(1, J, coeff_p);
 end
 
 if (priors.camera_prior == 1)
     % prior terms on rotations
    coeff_p = zeros((n_frames-1)*2, (K+6)*n_frames + K*3*n_points);
    for f=0:n_frames-1
        i_begin = f*(K+6)+K+1;
        j_begin = f*K+1;

        coeff_p(j_begin:j_begin+5, i_begin:i_begin+5) = 1;
        coeff_p(j_begin:j_begin+5, i_begin+(6+K):i_begin+(6+K)+5) = 1;

    end
    J = cat(1, J, coeff_p);
 end
 
J = sparse(J);

% % You can see easily the Jacobian pattern using the command spy(J)
disp('Observe the Jacobian pattern...')
spy(J)
