function [ FMatrix_deNormalized ] = deNorml( FMatrix, N_x, N_xp )

T=N_x;
U=N_xp;

FMatrix_deNormalized = T'*FMatrix*U;

end

