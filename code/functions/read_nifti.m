function [ x ] = read_nifti( y )

x = spm_read_vols(spm_vol(y));
end

